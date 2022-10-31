#include "global.h"
#include "Feistel.h"
#include "blake3.h"
#include <stdexcept>
#include <random>

Feistel::Feistel(int64_t domain_size, int num_rounds): domain_size{domain_size}, num_rounds{num_rounds}
{
    base = (int)ceil(log2(domain_size) / 2);
    std::random_device dev;
    seed = dev();
}

int64_t Feistel::encrypt(int64_t input)
{
    int64_t ander = (1 << base) - 1;
    int64_t l = input >> base;
    int64_t r = input & ander;
    int64_t tmp;
    int64_t buffer[2];
    blake3_hasher hasher;
    for (int i = 0; i < num_rounds; i++)
    {
        buffer[0] = r;
        buffer[1] = i ^ seed;
        blake3_hasher_init(&hasher);
        blake3_hasher_update(&hasher, buffer, 2 * sizeof(int64_t));
        blake3_hasher_finalize(&hasher, (uint8_t *)buffer, 2 * sizeof(int64_t));
        tmp = r;
        r = (l ^ buffer[1]) & ander;
        l = tmp;
    }
    return (l << base) | r;
}

int64_t Feistel::permute(int64_t input)
{
    if (input >= domain_size)
        throw std::invalid_argument("Index out of domain!");
    int64_t ret = encrypt(input);
    while (ret >= domain_size)
        ret = encrypt(ret);
    return ret;
}