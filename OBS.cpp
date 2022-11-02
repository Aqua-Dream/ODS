#include "OBS.h"
#include <cmath>
#include <stdexcept>
#include <random>
#include <omp.h>
#include <algorithm>
#include <boost/sort/sort.hpp>

// void RandomBinAssign(VectorSlice &vs, int64_t num_bins)
// {
//     // no need to add I/O cost, since it can actually be generated in memory at the first level in the fly
//     std::random_device dev;
//     std::vector<std::mt19937> rngs(NUM_THREADS);
//     std::uniform_int_distribution<int64_t> unif(0, num_bins - 1);
//     int64_t size = vs.size();
//     for (int i = 0; i < NUM_THREADS; i++)
//         rngs[i] = std::mt19937(dev());
// #pragma omp parallel for
//     for (int64_t i = 0; i < size; i++)
//         vs[i] = unif(rngs[omp_get_thread_num()]);
// }

ObliBucketSort::ObliBucketSort(IOManager &iom, int64_t dataSize, int blockSize, int sigma)
    : N{dataSize}, M{(int64_t)iom.GetInternalMemory().size()}, B{blockSize}, m_iom{iom}, m_intmem{iom.GetInternalMemory()}, sigma{sigma}
{
    if (N % B != 0 || M % B != 0)
        throw std::invalid_argument("N and M must be multiples of B!");
    double kappa = sigma * 0.693147;
    Z = 6 * (kappa + log(2.0 * N));
    Z = 6 * (kappa + log(2.0 * N / Z));
    Z = ceil(Z / B) * B;
}

void ObliBucketSort::RandomBinPermute(std::vector<int64_t> &input, std::vector<int64_t> &output)
{
    int n = ceil(log2(2.0 * N / Z));
    int m = n / 2;
    if (m > floor(log2(M / Z)))
        throw std::logic_error("Input size too large: not implemented yet!");
    FirstLevelPermute(input, output);
    SecondLevelPermute(output);
}

void ObliBucketSort::FirstLevelPermute(std::vector<int64_t> &input, std::vector<int64_t> &output)
{
    int n = ceil(log2(2.0 * N / Z));
    int m = n/ 2;
    int64_t num_bins = 1 << n;
    int64_t bins_perload = 1 << m;
    int64_t num_loads = 1 << (n - m);
    int64_t eles_perload = bins_perload * Z;
    int64_t real_eles_perload = ceil(N / num_bins) * bins_perload;
    output.resize(num_bins*Z);
    std::random_device dev;
    std::mt19937 rng(dev());
    for (int64_t i = 0; i < num_loads; i++)
    {
        VectorSlice dataext(input, i * real_eles_perload, real_eles_perload, true);
        int64_t num_reals = dataext.size();
        VectorSlice realint(m_intmem, 0, num_reals);
        VectorSlice intslice(m_intmem, 0, eles_perload);
        m_iom.DataTransfer(dataext, realint);
        std::fill(realint.end(), intslice.end(), DUMMY);
        std::shuffle(intslice.begin(), intslice.end(), rng);
        VectorSlice outslice(output, i * eles_perload, eles_perload);
        m_iom.DataTransfer(intslice, outslice);
    }
}

void ObliBucketSort::SecondLevelPermute(std::vector<int64_t> &output)
{
    int n = ceil(log2(2.0 * N / Z));
    int base = n/2;
    int m = n - base;
    int64_t num_bins = 1 << n;
    int64_t bins_perload = 1 << m;
    int64_t num_loads = 1 << base;
    int64_t eles_perload = bins_perload * Z;
    std::random_device dev;
    std::mt19937 rng(dev());
    for (int64_t i = 0; i < num_loads; i++)
    {
        for (int64_t j = 0; j < bins_perload; j++)
        {
            VectorSlice extslice_bin(output, Z * (i | (j << base)), Z);
            VectorSlice intslice_bin(m_intmem, j * Z, Z);
            m_iom.DataTransfer(extslice_bin, intslice_bin);
        }
        std::shuffle(m_intmem.begin(), m_intmem.begin() + eles_perload, rng);
        for (int64_t j = 0; j < bins_perload; j++)
        {
            VectorSlice intslice_bin(m_intmem, j * Z, Z);
            VectorSlice extslice_bin(output, Z * (i + (j << base)), Z);
            m_iom.DataTransfer(intslice_bin, extslice_bin);
        }
    }
}

void ObliBucketSort::Sort(std::vector<int64_t> &input, std::vector<int64_t> &output, bool printInfo)
{
    if (printInfo)
    {
        std::cout << "------------Bucket Sort ------------" << std::endl;
        Tick("Total");
        Tick("Random bin permutation");
    }
    RandomBinPermute(input, output);
    if (printInfo)
    {
        Tick("Random bin permutation");
        Tick("Final sort");
    }
    boost::sort::block_indirect_sort(output.begin(), output.end(), NUM_THREADS);
    m_iom.m_numIOs += 4*N/B; // merge sort, omitted
    output.resize(N);
    if (printInfo)
    {
        Tick("Final sort");
        Tick("Total");
        std::cout << "Num IOs: " << (float)m_iom.GetNumIOs() * B / N << "*N/B" << std::endl;
        m_iom.ClearIO();
    }
}