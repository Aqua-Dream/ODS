#include "memory.h"
#include "global.h"
#include <stdexcept>
#include <cassert>

void SmallSample(IOManager &iom, std::vector<uint64_t> &extmem, uint64_t blockSize, uint64_t sampleSize, bool oblivious)
{
    auto& intmem = iom.GetInternalMemory();
    uint64_t N = extmem.size();
    uint64_t M = intmem.size();
    uint64_t B = blockSize;
    //assert(N % B == 0);
    uint64_t n = sampleSize;
    if (n + B > M)
        throw std::invalid_argument("Sample does not fit into internal memory!");
    uint64_t num_sampled_eles = 0;
    for (uint64_t ifirst=0; ifirst<extmem.size(); ifirst+=B)
    {
        uint64_t num_sampled_eles_block = 0;
        for (uint64_t j = 0; j < B && N > 0; j++)
        {
            if (RandRange(0, N) < n)
            {
                intmem[num_sampled_eles + num_sampled_eles_block] = j;
                num_sampled_eles_block++;
                n--;
            }
            N--;
        }
        if (oblivious || num_sampled_eles_block > 0)
        {
            VectorSlice extslice(extmem, ifirst, B, true);
            VectorSlice intslice(intmem, M - B, B);
            iom.DataTransfer(extslice, intslice);
            for (uint64_t j = 0; j < num_sampled_eles_block; j++)
                intmem[num_sampled_eles + j] = intslice[intmem[num_sampled_eles + j]];
        }
        num_sampled_eles += num_sampled_eles_block;
    }
}