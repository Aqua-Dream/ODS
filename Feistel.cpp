#include "global.h"
#include "Feistel.h"
#include "blake3.h"
#include <stdexcept>

Feistel::Feistel(uint64_t domain_size, int num_rounds)
{
    this->domain_size = domain_size;
    this->base = (int)ceil(log2(domain_size) / 2);
    this->num_rounds = num_rounds;
}

uint64_t Feistel::encrypt(uint64_t input)
{
    uint64_t ander = (1 << base) - 1;
    uint64_t l = input >> base;
    uint64_t r = input & ander;
    uint64_t tmp;
    uint64_t buffer[2];
    blake3_hasher hasher;
    for (uint64_t i = 0; i < num_rounds; i++)
    {
        buffer[0] = r;
        buffer[1] = i;
        blake3_hasher_init(&hasher);
        blake3_hasher_update(&hasher, buffer, 2 * sizeof(uint64_t));
        blake3_hasher_finalize(&hasher, (uint8_t *)buffer, 2 * sizeof(uint64_t));
        tmp = r;
        r = l ^ buffer[1] & ander;
        l = tmp;
    }
    return (l << base) | r;
}

uint64_t Feistel::permute(uint64_t input)
{
    if (input >= domain_size)
        throw std::invalid_argument("Index out of domain!");
    uint64_t ret = encrypt(input);
    while (ret >= domain_size)
        ret = encrypt(ret);
    return ret;
}