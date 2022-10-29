#pragma once
#include "memory.h"
#include <vector>
#include<cstdint>

// the sample fits into internal memory
// oblivious: should the access pattern be independent of the sample ?
void SmallSample(IOManager &iom, std::vector<uint64_t> &externMem, uint64_t blockSize, uint64_t sampleSize, bool oblivious);

uint64_t InMemoryBernoulliSample(IOManager &iom, std::vector<uint64_t> &extmem, uint64_t blockSize, double sampleRate, bool oblivious);