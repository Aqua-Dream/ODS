#include "memory.h"
#include "global.h"
#include <stdexcept>
#include <cassert>
#include <omp.h>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <algorithm>

uint64_t InMemoryBernoulliSample(IOManager &iom, std::vector<uint64_t> &extmem, uint64_t blockSize, double sampleRate, bool oblivious)
{
    auto &intmem = iom.GetInternalMemory();
    uint64_t N = extmem.size();
    uint64_t M = intmem.size();
    uint64_t B = blockSize;
    double alpha = sampleRate;
    if (alpha * N > 0.9 * M)
        throw std::invalid_argument("Sample rate is too large: using in memory sample may fail.");
    uint64_t intPos = 0;
    boost::random_device dev;
    boost::random::binomial_distribution<int> binom(B, alpha);
    std::vector<boost::random::mt19937> rngs(NUM_THREADS);
    for (int i = 0; i < NUM_THREADS; i++)
        rngs[i] = boost::random::mt19937(dev());

    for (uint64_t extPos = 0; extPos < N; extPos += B * NUM_THREADS)
    {
        uint64_t num_block_remained = (N - extPos) / B;
        int num_threads = num_block_remained < NUM_THREADS ? num_block_remained : NUM_THREADS;
        if (intPos + num_threads * B > M)
        {
            std::cerr << "In memory sample overflow!" << std::endl;
            exit(-1);
        }
        uint64_t num_IOs = 0;

        #pragma omp parallel for reduction(+ : num_IOs) num_threads(num_threads)
        for (uint64_t j = 0; j < num_threads; j++)
        {
            int num_to_sample = binom(rngs[j]);
            VectorSlice intslice(intmem, M - (j + 1) * B, B);
            if (oblivious || num_to_sample > 0)
            {
                VectorSlice extslice(extmem, extPos + j * B, B);
                intslice.CopyDataFrom(extslice);
                std::random_shuffle(intslice.begin(), intslice.end());
            }
            std::fill(intslice.begin() + num_to_sample, intslice.end(), DUMMY);
            num_IOs += num_to_sample;
        }
        iom.m_numIOs += num_IOs;
        for (uint64_t j = M - 1; j >= M - num_threads * B; j--)
        {
            if (intmem[j] != DUMMY)
                intmem[intPos++] = intmem[j];
        }
    }
    return intPos;
}

void SmallSample(IOManager &iom, std::vector<uint64_t> &extmem, uint64_t blockSize, uint64_t sampleSize, bool oblivious)
{
    auto &intmem = iom.GetInternalMemory();
    uint64_t N = extmem.size();
    uint64_t M = intmem.size();
    uint64_t B = blockSize;
    if (M % B != 0 || N % B != 0)
        throw std::invalid_argument("N and M must be multiples of B!");
    uint64_t n = sampleSize;
    if (n + B > M)
        throw std::invalid_argument("Sample does not fit into internal memory!");
    uint64_t num_sampled_eles = 0;
    for (uint64_t ifirst = 0; ifirst < extmem.size(); ifirst += B)
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