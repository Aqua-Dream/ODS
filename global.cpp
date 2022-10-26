#include <random>
#include <cstdint>
#include <vector>
#include "global.h"

std::random_device dev;
std::mt19937 rng(0); 

// "end" not included
uint64_t RandRange(uint64_t start, uint64_t end)
{
    std::uniform_int_distribution<uint64_t>  distr(start, end-1);
    return distr(rng);
}
