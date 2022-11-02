#include "ODS.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include "Feistel.h"
#include <omp.h>
#include <iostream>
#include <boost/sort/sort.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/pow.hpp>

struct equation_to_solve
{
    // Functor returning both 1st and 2nd derivatives.
    equation_to_solve(double a) : m_arg(a) {}

    std::tuple<double, double, double> operator()(double x)
    {
        using boost::math::pow;
        // Return both f(x) and f'(x) and f''(x).
        double fx = pow<2>(x) * (1 - x) - pow<2>(1 + 2 * x) * (1 + x / 2) * pow<2>(1 + x) * m_arg;
        double dx = (2 - 3 * x) * x - m_arg * (6.5 + x * (34 + x * (55.5 + x * (40 + x * 10))));
        double d2x = 2 - 6 * x - m_arg * (34 + x * (111 + x * 120 + x * 40));
        return std::make_tuple(fx, dx, d2x); // 'return' fx, dx and d2x.
    }

private:
    double m_arg; // to be 'fifth_rooted'.
};

// return -1 if no solution exists
// For one-level, a=2NB/M^2 *(kappa+1+2*log(N/M))
// For two-level, a=2B*sqrt(N/M)/M *(kappa+2+1.5*log(N/M))
double argsolver(double a)
{
    // find beta
    using namespace std;                // Help ADL of std functions.
    using namespace boost::math::tools; // for halley_iterate.
    double guess = 0.2;                 // Rough guess is to divide the exponent by five.
    double min = 0.0;                   // Minimum possible value is half our guess.
    double max = 1.0;                   // Maximum possible value is twice our guess.
    const int digits = 5;
    const boost::uintmax_t maxit = 10;
    boost::uintmax_t it = maxit;
    double result = halley_iterate(equation_to_solve(a), guess, min, max, digits, it);
    if (result < 1e-5 || result > 1 - 1e-5) // fail to solve
        result = -1;
    return result;
}

// return the q-quantile of the first n elements in data
// length of quantile: q-1
std::vector<int64_t> GetQuantile(VectorSlice vs, int q)
{
    int64_t n = vs.size();
    int64_t step = n / q;
    int remain = n % q;
    std::vector<int64_t> pivots(q - 1);
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
    int64_t i = 0;
    int64_t j = data.size() - 1;
    while (1)
    {
        while (i < j && data[i] != DUMMY)
            i++; // find the first dummy from left
        while (i < j && data[j] == DUMMY)
            j--; // find the first non-dummy from right
        if (i >= j)
            break;
        int64_t tmp = data[i];
        data[i] = data[j];
        data[j] = tmp;
    }
    return VectorSlice(data, 0, i);
}

ObliDistSort::ObliDistSort(IOManager &iom, int64_t dataSize, int blockSize, int sigma)
    : N{dataSize}, M{(int64_t)iom.GetInternalMemory().size()}, B{blockSize}, m_iom{iom}, m_intmem{iom.GetInternalMemory()}, sigma{sigma}
{
    if (M % B != 0 || N % B != 0)
        throw std::invalid_argument("N and M must be multiples of B!");
    alpha = beta = p0 = p = -1; // not implemented
}

int64_t ObliDistSort::GetSampleSizeEachMemload()
{
    return 1.2 * alpha * M; 
}

int64_t ObliDistSort::GetSampleSize()
{
    return GetSampleSizeEachMemload() * ceil((float)N / M); 
}

// each element is sampled with probability alpha, independently
void ObliDistSort::Sample(std::vector<int64_t> &input, std::vector<int64_t> &output, SortType sorttype)
{
    int64_t outPos = 0;
    output.resize(GetSampleSize());
    std::random_device dev;
    std::binomial_distribution<int> binom(B, alpha);
    std::vector<std::mt19937> rngs(NUM_THREADS);
    for (int i = 0; i < NUM_THREADS; i++)
        rngs[i] = std::mt19937(dev());
    for (int64_t extPos = 0; extPos < N; extPos += M)
    {
        int64_t memload = N - extPos;
        if (memload > M)
            memload = M;
        int64_t num_IOs = 0;
        std::fill(m_intmem.begin(), m_intmem.end(), DUMMY);
#pragma omp parallel for reduction(+ \
                                   : num_IOs)
        for (int64_t j = 0; j < memload / B; j++)
        {
            std::mt19937 &rng = rngs[omp_get_thread_num()];
            int num_to_sample = binom(rng);
            VectorSlice intslice(m_intmem, j * B, B);
            if (sorttype == SortType::TIGHT || num_to_sample > 0)
            {
                num_IOs++;
                VectorSlice extslice(input, extPos + j * B, B);
                intslice.CopyDataFrom(extslice);
                // randomly take "num_to_sample" elements from the block
                for (int k = 0; k < num_to_sample; k++)
                {
                    int64_t tmp = intslice[k];
                    std::uniform_int_distribution<int> unif(k, num_to_sample - 1);
                    int l = unif(rng);
                    intslice[k] = intslice[l];
                    intslice[l] = tmp;
                }
                std::fill(intslice.begin() + num_to_sample, intslice.end(), DUMMY);
            }
        }
        m_iom.m_numIOs += num_IOs;
        VectorSlice intslice(m_intmem, 0, memload);
        auto realslice = Compact(intslice);
        int64_t outsize = realslice.size();
        if (sorttype == SortType::TIGHT)
        {
            outsize = GetSampleSizeEachMemload();
            if (realslice.size() > outsize)
            {
                std::cerr << "Data overflow when sampling!" << std::endl;
                exit(-1);
            }
        }
        VectorSlice outslice(output, outPos, outsize);
        m_iom.DataTransfer(realslice, outslice);
        outPos += outsize;
    }
    int64_t dummyToPad = outPos % B > 0 ? B - outPos % B : 0;
    output.resize(outPos + dummyToPad);
    std::fill_n(output.begin() + outPos, dummyToPad, DUMMY);
}

std::vector<int64_t> ObliDistSort::GetPivots(std::vector<int64_t> &data, SortType sorttype)
{
    throw std::logic_error("Sorting for general ODS has not been implemented yet!");
}

std::vector<std::vector<int64_t> *> ObliDistSort::Partition(std::vector<int64_t> &data, std::vector<int64_t> &pivots, bool isFirstLevel)
{
    int64_t memload = isFirstLevel ? (ceil(M / (1 + 2 * beta) / B) * B) : M; // be the multiple of B
    int64_t data_size = data.size();
    if (data_size % B != 0)
        throw std::logic_error("Data size must be the multiple of block size!");
    int64_t num_memloads = ceil((float)data_size / memload);
    int num_buckets = pivots.size() + 1;
    int64_t unit = ceil(float(M) / num_buckets);
    std::vector<uint64_t> posList(num_buckets - 1);
    std::vector<std::vector<int64_t> *> buckets(num_buckets);
    int64_t bucket_size = ceil((float)(unit * num_memloads) / B) * B;
    for (int i = 0; i < num_buckets; i++)
        buckets[i] = new std::vector<int64_t>(bucket_size, DUMMY);
    Feistel fs(data_size / B);
    for (int64_t i = 0; i < num_memloads; i++)
    {
        int64_t actual_load = (data_size < (i + 1) * memload) ? (data_size - i * memload) : memload;
        int64_t blocks_thisload = actual_load / B;
        VectorSlice intslice(m_intmem, 0, actual_load);
#pragma omp parallel for
        for (int64_t j = 0; j < blocks_thisload; j++)
        {
            int64_t blockId = j + i * memload / B;
            if (isFirstLevel)
                blockId = fs.permute(blockId);
            std::copy_n(data.begin() + blockId * B, B, m_intmem.begin() + j * B);
        }
        m_iom.m_numIOs += blocks_thisload;
        boost::sort::block_indirect_sort(intslice.begin(), intslice.end(), NUM_THREADS);
        // remove dummy
        VectorSlice _intslice(intslice, intslice.begin(), std::lower_bound(intslice.begin(), intslice.end(), DUMMY));
        std::vector<int64_t>::iterator iprevPos = _intslice.begin();
#pragma omp parallel
        {
#pragma omp for
            for (int j = 0; j < num_buckets - 1; j++)
            {
                auto pivot = pivots[j];
                posList[j] = std::distance(_intslice.begin(), std::lower_bound(_intslice.begin(), _intslice.end(), pivot));
            }
#pragma omp for
            for (int j = 0; j < num_buckets; j++)
            {
                int64_t prevPos = (j == 0) ? 0 : posList[j - 1];
                int64_t endPos = (j == num_buckets - 1) ? _intslice.size() : posList[j];
                VectorSlice intBucket(_intslice, prevPos, endPos - prevPos);
                VectorSlice extBucket(*buckets[j], unit * i, unit);
                if (intBucket.size() > unit)
                {
                    std::cerr << "Bucket overflow!" << std::endl;
                    exit(-1);
                }
                extBucket.CopyDataFrom(intBucket, true);
            }
        }
        m_iom.m_numIOs += num_buckets * ceil((float)unit / B);
    }
    return buckets;
}

void ObliDistSort::FinalSorting(std::vector<std::vector<int64_t> *> &buckets, std::vector<int64_t> &out, SortType sorttype)
{
    int64_t bucket_size = buckets[0]->size();
    if (bucket_size > M)
        throw std::invalid_argument("Memory not enough for final sorting!");
    int64_t num_buckets = buckets.size();
    if (sorttype == TIGHT)
        out.resize(N);
    else
        out.resize(bucket_size * num_buckets);
    int64_t pos = 0;
    for (int64_t i = 0; i < num_buckets; i++)
    {
        VectorSlice extslice(*buckets[i], 0, bucket_size);
        VectorSlice intslice(m_intmem, 0, bucket_size);
        m_iom.DataTransfer(extslice, intslice);
        boost::sort::block_indirect_sort(intslice.begin(), intslice.end(), NUM_THREADS);
        VectorSlice _intslice(intslice, intslice.begin(), std::lower_bound(intslice.begin(), intslice.end(), DUMMY));
        int64_t num_to_write = (sorttype == TIGHT) ? _intslice.size() : bucket_size;
        VectorSlice outslice(out, pos, num_to_write);
        m_iom.DataTransfer(_intslice, outslice);
        pos += num_to_write;
    }
    out.resize(pos);
}

void ObliDistSort::Sort(std::vector<int64_t> &input, std::vector<int64_t> &output, SortType sorttype)
{
    throw std::logic_error("Sorting for general ODS has not been implemented yet!");
}
