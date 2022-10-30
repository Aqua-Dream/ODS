#pragma once
#include <random>
#include <cstdint>
#include <cmath>


class Feistel
{
public:
    Feistel(int64_t domain_size, int num_rounds = 3);
    int64_t permute(int64_t input);

private:
    int64_t domain_size;
    int base, num_rounds;
    int64_t encrypt(int64_t input);
};