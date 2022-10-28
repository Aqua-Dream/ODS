#pragma once

#include <cstdint>
#include "memory.h"
#include "global.h"

//void Quantile(std::vector<uint64_t>& data, uint64_t n, int q); // for test only

//void InternalPartition(VectorSlice &data, VectorSlice &pivots, VectorSlice &posList);


class OneLevel
{
public:
    enum SortType{ TIGHT,LOOSE };
    // failure probability bounded by 2^(-sigma)
    OneLevel(IOManager &iom, uint64_t dataSize, uint64_t blockSize, int sigma);
    void GetPivots(std::vector<uint64_t> &extint, SortType sorttype);
    // extout should not be the same with extin
    std::vector<std::vector<uint64_t>*> FirstLevelPartition(std::vector<uint64_t> &extint);
    void FinalSorting(std::vector<std::vector<uint64_t>*> &buckets, std::vector<uint64_t> &out, SortType sorttype);
    void Sort(std::vector<uint64_t> &extint, std::vector<uint64_t> &extout, SortType sorttype);
private:
    uint64_t N, M, B;
    double alpha, beta;
    int p;
    IOManager &m_iom;
    std::vector<uint64_t>& m_intmem;
    
};