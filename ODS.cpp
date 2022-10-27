#include "ODS.h"
#include <cmath>
#include <stdexcept>
#include "sampling.h"
#include "argsolver.h"
#include <algorithm>
#include <cassert>
#include "Feistel.h"

// compute the q-quantile of the first n elements in data, and put them at the end of data
// length of quantile: q-1
void Quantile(std::vector<uint64_t> &data, uint64_t n, int q)
{
    uint64_t M = data.size();
    uint64_t step = n / q;
    int remain = n % q;
    uint64_t pivotId = M - 1;
    // dataIdx<n<=M<int_max, so there is no overflow
    int dataIdx = n;
    int i = 0;
    for (; i < remain; i++)
    {
        dataIdx -= step + 1;
        data[pivotId--] = data[dataIdx];
    }
    for (; i < q - 1; i++)
    {
        dataIdx -= step;
        data[pivotId--] = data[dataIdx];
    }
    assert(dataIdx == step);
    assert(pivotId == M - q);
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

// return: the number of elements smaller than pivot
uint64_t InternalPartitionLevel(VectorSlice &data, uint64_t pivot)
{
    uint64_t i = 0;
    for (uint64_t j = 0; j < data.size(); ++j)
    {
        if (pivot > data[j])
        {
            uint64_t tmp = data[i];
            data[i] = data[j];
            data[j] = tmp;
            i++;
        }
    }
    return i;
}

void InternalPartition(VectorSlice &data, VectorSlice &pivots, VectorSlice &posList)
{
    uint64_t pivotIdx, pivot, dataIdx;
    pivotIdx = pivots.size() >> 1;
    pivot = pivots[pivotIdx];
    dataIdx = InternalPartitionLevel(data, pivot);
    posList[pivotIdx] = data.GetStartPos() + dataIdx;
    if (dataIdx > 0 && pivotIdx > 0)
    {
        VectorSlice dataLeft(data, 0, dataIdx);
        VectorSlice pivotsLeft(pivots, 0, pivotIdx);
        VectorSlice posLeft(posList, 0, pivotIdx);
        InternalPartition(dataLeft, pivotsLeft, posLeft);
    }
    else if (dataIdx <= 0)
    {
        for (uint64_t i = 0; i < pivotIdx; i++)
            posList[i] = data.GetStartPos();
    }
    if (dataIdx < data.size() && pivotIdx + 1 < pivots.size())
    {
        VectorSlice dataRight(data, dataIdx, data.size() - dataIdx);
        VectorSlice pivotsRight(pivots, pivotIdx + 1, pivots.size() - pivotIdx - 1);
        VectorSlice posRight(posList, pivotIdx + 1, pivots.size() - pivotIdx - 1);
        InternalPartition(dataRight, pivotsRight, posRight);
    }
    else if (dataIdx >= data.size())
    {
        for (uint64_t i = pivotIdx + 1; i < posList.size(); i++)
            posList[i] = data.GetEndPos();
    }
}

OneLevel::OneLevel(IOManager &iom, uint64_t dataSize, uint64_t blockSize, int sigma)
    : N{dataSize}, M{iom.GetInternalMemory().size()}, B{blockSize}, m_iom{iom}, m_intmem{iom.GetInternalMemory()}
{
    if (M % B != 0 || N % B != 0)
        throw std::invalid_argument("N and M must be multiple of B!");
    double kappa = sigma * 0.693147;
    // has a solution iff a<0.012858
    double a = 2.0 * N * B / M * (kappa + 1 + 2 * log(N / M)) / M;
    beta = argsolver(a);
    if (beta == -1)
        throw std::invalid_argument("Invalid parameters for one-level sorting!");
    alpha = (kappa + 1 + log(N)) * 4 * (1 + 1 / beta) * (1 / beta + 2) / M;
    while (alpha * N >= M - B)
    {
        beta *= 2;
        alpha = (kappa + 1 + log(N)) * 4 * (1 + 1 / beta) * (1 / beta + 2) / M;
    }
    if (beta >= 1)
        throw std::invalid_argument("Invalid parameters for one-level sorting!");
    p = (int)ceil((1 + 2 * beta) * N / M);
}

void OneLevel::DecidePivots(std::vector<uint64_t> &extint, SortType sorttype)
{
    uint64_t n = (int)ceil(alpha * N);
    SmallSample(m_iom, extint, B, n, sorttype == TIGHT);
    std::sort(m_intmem.begin(), m_intmem.begin() + n);
    // move pivots to the end of the memory
    Quantile(m_intmem, n, p);
}

std::vector<std::vector<uint64_t> *> OneLevel::FirstLevelPartition(std::vector<uint64_t> &extint)
{
    uint64_t memload = ceil(M / (1 + 2 * beta) / B) * B; // be the multiple of B
    uint64_t num_memloads = ceil((float)N / memload);
    uint64_t unit = ceil(float(M) / p);
    if (memload + 2 * (p - 1) > M)
        throw std::invalid_argument("Memory not enough!");
    std::vector<std::vector<uint64_t> *> buckets(p);
    for (int i = 0; i < p; i++)
        buckets[i] = new std::vector<uint64_t>(unit * num_memloads);
    VectorSlice pivots(m_intmem, M - p + 1, p - 1);
    VectorSlice posList(m_intmem, M - 2 * (p - 1), p - 1);
    Feistel fs(N / B);
    for (uint64_t i = 0; i < num_memloads; i++)
    {
        uint64_t actual_load = (N < (i + 1) * memload) ? (N - i * memload) : memload;
        for (uint64_t j = 0; j < actual_load; j++)
        {
            uint64_t blockId = fs.permute(j + i * memload);
            VectorSlice extslice(extint, blockId * B, B);
            VectorSlice intslice(m_intmem, j * B, B);
            m_iom.DataTransfer(extslice, intslice);
        }
        VectorSlice intdata(m_intmem, 0, actual_load);
        InternalPartition(intdata, pivots, posList);
        for (uint64_t j = 0; j < p; j++)
        {
            uint64_t prevPos = (j == 0) ? 0 : posList[j - 1];
            uint64_t endPos = (j == p - 1) ? actual_load : posList[j];
            VectorSlice intBucket(m_intmem, prevPos, endPos - prevPos);
            VectorSlice extBucket(*buckets[j], unit * i, unit);
            m_iom.DataTransfer(intBucket, extBucket);
        }
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
        std::sort(realslice.begin(), realslice.end());
        uint64_t num_to_write = (sorttype == TIGHT) ? realslice.size() : bucket_size;
        VectorSlice outslice(out, pos, num_to_write);
        m_iom.DataTransfer(realslice, outslice);
        pos += num_to_write;
    }
}

void OneLevel::Sort(std::vector<uint64_t> &extint, std::vector<uint64_t> &extout, SortType sorttype)
{
    DecidePivots(extint, sorttype);
    auto buckets = FirstLevelPartition(extint);
    FinalSorting(buckets, extout, sorttype);
    for(auto vs: buckets)
        delete vs;
}