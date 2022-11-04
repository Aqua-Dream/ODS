#pragma once

#include <cstdint>
#include "memory.h"
#include "global.h"

double argsolver(double a);

// return the q-quantile of the first n elements in data
// length of quantile: q-1
template <typename T>
std::vector<T> GetQuantile(VectorSlice<T> vs, int q)
{
    int64_t n = vs.size();
    int64_t step = n / q;
    int remain = n % q;
    std::vector<T> pivots(q - 1);
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

enum SortType
{
    TIGHT,
    LOOSE
};

template <typename T>
class ObliDistSort
{
public:
    // failure probability bounded by 2^(-sigma)
    ObliDistSort(IOManager<T> &iom, int64_t dataSize, int blockSize, int sigma);
    void Sort(std::vector<T> &input, std::vector<T> &output, SortType sorttype);

protected:
    int64_t N, M;
    int B;
    double alpha, beta;
    int p0, p, sigma, num_levels;
    IOManager<T> &m_iom;
    std::vector<T> &m_intmem;
    int64_t GetSampleSizeEachMemload();
    int64_t GetSampleSize();
    void Sample(std::vector<T> &input, std::vector<T> &output, SortType sorttype);
    std::vector<T> GetPivots(std::vector<T> &data, SortType sorttype);
    std::vector<std::vector<T> *> Partition(std::vector<T> &data, std::vector<T> &pivots, bool isFirstLevel);
    void FinalSorting(std::vector<std::vector<T> *> &buckets, std::vector<T> &out, SortType sorttype);
};

template <typename T>
class OneLevel : public ObliDistSort<T>
{

public:
    OneLevel(IOManager<T> &iom, int64_t dataSize, int blockSize, int sigma);
    void Sort(std::vector<T> &input, std::vector<T> &output, SortType sorttype, bool printInfo=false);
private:
    std::vector<T> GetPivots(std::vector<T> &data, SortType sorttype);
};

template <typename T>
class TwoLevel : public ObliDistSort<T>
{
public:
    TwoLevel(IOManager<T> &iom, int64_t dataSize, int blockSize, int sigma);
    void Sort(std::vector<T> &input, std::vector<T> &output, SortType sorttype, bool printInfo=false);
private:
    std::vector<T> GetPivots(std::vector<T> &data, SortType sorttype);
};