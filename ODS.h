#pragma once

#include <cstdint>
#include "memory.h"
#include "global.h"

void Quantile(std::vector<uint64_t>& data, uint64_t n, int q); // for test only

void InternalPartition(VectorSlice &data, VectorSlice &pivots, VectorSlice &posList);


class OneLevel
{
public:
    enum Type{ TIGHT,LOOSE };
    // failure probability bounded by 2^(-sigma)
    OneLevel(IOManager &iom, std::vector<uint64_t> & extdata, uint64_t blockSize, int sigma);
    void DecidePivots(Type m_type);
    void FirstLevelPartition();
    
private:
    uint64_t N, M, B;
    double alpha, beta;
    int p;
    IOManager &m_iom;
    std::vector<uint64_t>& m_intmem;
    std::vector<uint64_t>& m_extdata;
    
};