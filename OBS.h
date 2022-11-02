#pragma once

#include <cstdint>
#include "memory.h"
#include "global.h"

class ObliBucketSort
{
public:
    // failure probability bounded by 2^(-sigma)
    ObliBucketSort(IOManager &iom, int64_t dataSize, int blockSize, int sigma);
    void Sort(std::vector<int64_t> &input, std::vector<int64_t> &output, bool printInfo);

protected:
    int64_t N, M;
    int B, Z;
    int sigma, num_levels;
    IOManager &m_iom;
    std::vector<int64_t> &m_intmem;
    void FirstLevelPermute(std::vector<int64_t> &input, std::vector<int64_t> &output);
    void SecondLevelPermute(std::vector<int64_t> &output);
};
