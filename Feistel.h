#pragma once
#include <random>
#include <cstdint>
#include <cmath>


class Feistel
{
public:
    Feistel(uint64_t domain_size, int num_rounds = 3);
    uint64_t encrypt(uint64_t input);

private:
    uint64_t domain_size;
    int base, num_rounds;
};