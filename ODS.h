#pragma once

#include <cstdint>
#include "memory.h"
#include "global.h"

double argsolver(double a);
std::vector<int64_t> GetQuantile(VectorSlice vs, int q);

class ObliDistSort
{
public:
    enum SortType
    {
        TIGHT,
        LOOSE
    };
    // failure probability bounded by 2^(-sigma)
    ObliDistSort(IOManager &iom, int64_t dataSize, int blockSize, int sigma);
    void Sort(std::vector<int64_t> &input, std::vector<int64_t> &output, SortType sorttype);

protected:
    int64_t N, M;
    int B;
    double alpha, beta;
    int p0, p, sigma, num_levels;
    IOManager &m_iom;
    std::vector<int64_t> &m_intmem;
    void Sample(std::vector<int64_t> &input, std::vector<int64_t> &output, SortType sorttype);
    std::vector<int64_t> GetPivots(std::vector<int64_t> &data, SortType sorttype);
    std::vector<std::vector<int64_t> *> Partition(std::vector<int64_t> &data, std::vector<int64_t> &pivots, bool isFirstLevel);
    void FinalSorting(std::vector<std::vector<int64_t> *> &buckets, std::vector<int64_t> &out, SortType sorttype);
};

class OneLevel : public ObliDistSort
{
public:
    OneLevel(IOManager &iom, int64_t dataSize, int blockSize, int sigma);
    void Sort(std::vector<int64_t> &input, std::vector<int64_t> &output, SortType sorttype);
private:
    std::vector<int64_t> GetPivots(std::vector<int64_t> &data, SortType sorttype);
};

class TwoLevel : public ObliDistSort
{
public:
    TwoLevel(IOManager &iom, int64_t dataSize, int blockSize, int sigma);
    void Sort(std::vector<int64_t> &input, std::vector<int64_t> &output, SortType sorttype);
private:
    std::vector<int64_t> GetPivots(std::vector<int64_t> &data, SortType sorttype);
};