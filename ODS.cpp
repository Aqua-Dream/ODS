#include "ODS.h"
#include <cmath>
#include <stdexcept>
#include "sampling.h"
#include "argsolver.h"
#include <algorithm>
#include <cassert>
#include "Feistel.h"
#include <omp.h>
#include <iostream>
#include <boost/sort/sort.hpp>

// return the q-quantile of the first n elements in data
// length of quantile: q-1
std::vector<uint64_t> GetQuantile(VectorSlice vs, int q)
{
    uint64_t n = vs.size();
    uint64_t step = n / q;
    int remain = n % q;
    std::vector<uint64_t> pivots(q - 1);
    int dataIdx = -1;
    int i = 0;
    for (; i < remain; i++)
    {
        dataIdx += step + 1;
        pivots[i] = vs[dataIdx];
    }
    for (; i < q - 1; i++)
    {
        dataIdx += step;
        pivots[i] = vs[dataIdx];
    }
    return pivots;
}

// move dummy elements to the end; return the slice of real elements
VectorSlice Compact(VectorSlice data)
{
    uint64_t i = 0;
    uint64_t j = data.size() - 1;
    while (1)
    {
        while (i < j && data[i] != DUMMY)
            i++; // find the first dummy from left
        while (i < j && data[j] == DUMMY)
            j--; // find the first non-dummy from right
        if (i >= j)
            break;
        uint64_t tmp = data[i];
        data[i] = data[j];
        data[j] = tmp;
    }
    return VectorSlice(data, 0, i);
}

OneLevel::OneLevel(IOManager &iom, uint64_t dataSize, uint64_t blockSize, int sigma)
    : N{dataSize}, M{iom.GetInternalMemory().size()}, B{blockSize}, m_iom{iom}, m_intmem{iom.GetInternalMemory()}
{
    if (M % B != 0 || N % B != 0)
        throw std::invalid_argument("N and M must be multiples of B!");
    double kappa = sigma * 0.693147;
    // has a solution iff a<0.012858
    double a = 2.0 * N * B / M * (kappa + 1 + 2 * log(N / M)) / M;
    beta = argsolver(a);
    if (beta == -1)
        throw std::invalid_argument("Invalid parameters for one-level sorting!");
    alpha = (kappa + 1 + log(N)) * 4 * (1 + 1 / beta) * (1 / beta + 2) / M;
    while (alpha >= 1)
    {
        beta *= 2;
        alpha = (kappa + 1 + log(N)) * 4 * (1 + 1 / beta) * (1 / beta + 2) / M;
    }
    if (beta >= 1)
        throw std::invalid_argument("Invalid parameters for one-level sorting!");
    p = (int)ceil((1 + 2 * beta) * N / M);
    uint64_t memload = ceil(M / (1 + 2 * beta) / B) * B; // be the multiple of B
    uint64_t num_memloads = ceil((float)N / memload);
    uint64_t unit = ceil(float(M) / p);
    if (unit * num_memloads > M)
        p++;
}

void OneLevel::Sample(std::vector<uint64_t> &extin, std::vector<uint64_t> &extout, SortType sorttype)
{
    uint64_t outPos = 0;
    extout.resize((uint64_t)(1.5 * alpha * N));
    std::random_device dev;
    std::binomial_distribution<int> binom(B, alpha);
    std::vector<std::mt19937> rngs(NUM_THREADS);
    for (int i = 0; i < NUM_THREADS; i++)
        rngs[i] = std::mt19937(dev());
    for (uint64_t extPos = 0; extPos < N; extPos += M)
    {
        uint64_t memload = N - extPos;
        if (memload > M)
            memload = M;
        uint64_t num_IOs = 0;
        std::fill(m_intmem.begin(), m_intmem.end(), DUMMY);
#pragma omp parallel for reduction(+ \
                                   : num_IOs) num_threads(NUM_THREADS)
        for (uint64_t j = 0; j < memload / B; j++)
        {
            uint64_t num_to_sample = binom(rngs[omp_get_thread_num()]);
            VectorSlice intslice(m_intmem, j * B, B);
            if (sorttype == SortType::TIGHT || num_to_sample > 0)
            {
                num_IOs++;
                VectorSlice extslice(extin, extPos + j * B, B);
                intslice.CopyDataFrom(extslice);
                std::random_shuffle(intslice.begin(), intslice.end());
                std::fill(intslice.begin() + num_to_sample, intslice.end(), DUMMY);
            }
        }
        m_iom.m_numIOs += num_IOs;
        VectorSlice intslice(m_intmem, 0, memload);
        auto realslice = Compact(intslice);
        uint64_t outsize =  realslice.size();
        if (sorttype == SortType::TIGHT)
        {
            outsize = 1.5 * alpha * M;
            if (realslice.size() > outsize)
            {
                std::cerr << "Data overflow when sampling!" << std::endl;
                exit(-1);
            }

        }
        VectorSlice outslice (extout, outPos, outsize);
        m_iom.DataTransfer(realslice, outslice);
        outPos += outsize;
    }

    extout.resize(outPos);
}

std::vector<uint64_t> OneLevel::GetPivots(std::vector<uint64_t> &extint, SortType sorttype)
{
    std::vector<uint64_t> sample;
    Sample(extint, sample, sorttype);
    if(sample.size()>M)
        throw std::invalid_argument("Sample does not fit into memory!");
    VectorSlice sampleslice (sample, 0, sample.size());
    VectorSlice intslice(m_intmem, 0, sample.size());
    m_iom.DataTransfer(sampleslice, intslice);
    boost::sort::block_indirect_sort(intslice.begin(), intslice.end(), NUM_THREADS);
    auto realnum = std::distance(intslice.begin(), std::lower_bound(intslice.begin(),intslice.end(), DUMMY-1));
    // move pivots to the end of the memory
    VectorSlice sampleforquantile(intslice, 0, realnum);
    return GetQuantile(sampleforquantile, p);
}

std::vector<std::vector<uint64_t> *> OneLevel::FirstLevelPartition(std::vector<uint64_t> &extint, std::vector<uint64_t> &pivots)
{
    uint64_t memload = ceil(M / (1 + 2 * beta) / B) * B; // be the multiple of B
    uint64_t num_memloads = ceil((float)N / memload);
    uint64_t unit = ceil(float(M) / p);
    if (memload + 2 * (p - 1) > M)
        throw std::invalid_argument("Memory not enough!");
    std::vector<std::vector<uint64_t> *> buckets(p);
    for (int i = 0; i < p; i++)
        buckets[i] = new std::vector<uint64_t>(unit * num_memloads);
    Feistel fs(N / B);
    std::vector<uint64_t> posList(p - 1);
    for (uint64_t i = 0; i < num_memloads; i++)
    {
        uint64_t actual_load = (N < (i + 1) * memload) ? (N - i * memload) : memload;
        uint64_t blocks_thisload = actual_load / B;
#pragma omp parallel for num_threads(NUM_THREADS)
        for (uint64_t j = 0; j < blocks_thisload; j++)
        {
            uint64_t blockId = fs.permute(j + i * memload / B);
            std::copy_n(extint.begin() + blockId * B, B, m_intmem.begin() + j * B);
        }
        m_iom.m_numIOs += blocks_thisload;
        VectorSlice intdata(m_intmem, 0, actual_load);
        boost::sort::block_indirect_sort(intdata.begin(), intdata.end(), NUM_THREADS);
#pragma omp parallel for num_threads(NUM_THREADS)
        for (int j = 0; j < p - 1; j++)
            posList[j] = std::distance(intdata.begin(), std::lower_bound(intdata.begin(), intdata.end(), pivots[j]));
#pragma omp parallel for num_threads(NUM_THREADS)
        for (int j = 0; j < p; j++)
        {
            uint64_t prevPos = (j == 0) ? 0 : posList[j - 1];
            uint64_t endPos = (j == p - 1) ? intdata.size() : posList[j];
            VectorSlice intBucket(intdata, prevPos, endPos - prevPos);
            VectorSlice extBucket(*buckets[j], unit * i, unit);
            if (intBucket.size() > unit)
            {
                std::cerr << "Bucket overflow!" << std::endl;
                exit(-1);
            }
            extBucket.CopyDataFrom(intBucket, true);
        }
        m_iom.m_numIOs += ceil((float)p * unit / B);
    }
    return buckets;
}

void OneLevel::FinalSorting(std::vector<std::vector<uint64_t> *> &buckets, std::vector<uint64_t> &out, SortType sorttype)
{
    uint64_t bucket_size = buckets[0]->size();
    if (bucket_size > M)
        throw std::invalid_argument("Memory not enough for final sorting!");
    uint64_t num_buckets = buckets.size();
    if (sorttype == TIGHT)
        out.resize(N);
    else
        out.resize(bucket_size * num_buckets);
    uint64_t pos = 0;
    for (uint64_t i = 0; i < num_buckets; i++)
    {
        VectorSlice extslice(*buckets[i], 0, bucket_size);
        VectorSlice intslice(m_intmem, 0, bucket_size);
        m_iom.DataTransfer(extslice, intslice);
        auto realslice = Compact(intslice);
        boost::sort::block_indirect_sort(realslice.begin(), realslice.end(), NUM_THREADS);
        uint64_t num_to_write = (sorttype == TIGHT) ? realslice.size() : bucket_size;
        VectorSlice outslice(out, pos, num_to_write);
        m_iom.DataTransfer(realslice, outslice);
        pos += num_to_write;
    }
}

void OneLevel::Sort(std::vector<uint64_t> &extint, std::vector<uint64_t> &extout, SortType sorttype)
{
    Tick("Decide pivots");
    auto pivots = GetPivots(extint, sorttype);
    Tick("Decide pivots");

    Tick("Partition");
    auto buckets = FirstLevelPartition(extint, pivots);
    Tick("Partition");

    Tick("Final sorting");
    FinalSorting(buckets, extout, sorttype);
    Tick("Final sorting");
    for (auto vs : buckets)
        delete vs;
}