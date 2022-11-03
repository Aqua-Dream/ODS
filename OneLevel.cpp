#include "ODS.h"
#include <cmath>
#include <stdexcept>
#include <boost/sort/sort.hpp>

OneLevel::OneLevel(IOManager &iom, int64_t dataSize, int blockSize, int sigma)
    : ObliDistSort(iom, dataSize, blockSize, sigma)
{
    num_levels = 1;
    double kappa = sigma * 0.693147;
    // has a solution iff a<0.012858
    double a = 2.0 * N * B / M * (kappa + 1 + 2 * log(N / M)) / M;
    beta = argsolver(a);
    if (beta <= 0)
        throw std::invalid_argument("Invalid parameters for one-level sorting!");
    alpha = (kappa + 1 + log(N)) * 4 * (1 + 1 / beta) * (1 / beta + 2) / M;
    while (GetSampleSize() > M || alpha >= 2*beta)
    {
        beta *= 2;
        alpha = (kappa + 1 + log(N)) * 4 * (1 + 1 / beta) * (1 / beta + 2) / M;
    }
    if (beta >= 1)
        throw std::invalid_argument("Invalid parameters for one-level sorting!");
    p0 = (int)ceil((1 + 2 * beta) * N / M);
    int64_t memload = ceil(M / (1 + 2 * beta) / B) * B; // be the multiple of B
    int64_t num_memloads = ceil_divide(N, memload);
    int64_t unit = ceil(float(M) / p0);
    if (unit * num_memloads > M)
        p0++;
}

std::vector<int64_t> OneLevel::GetPivots(std::vector<int64_t> &data, SortType sorttype)
{
    std::vector<int64_t> sample;
    Sample(data, sample, sorttype);
    if (sample.size() > M)
        throw std::invalid_argument("Sample does not fit into memory!");
    VectorSlice sampleslice(sample, 0, sample.size());
    VectorSlice intslice(m_intmem, 0, sample.size());
    m_iom.DataTransfer(sampleslice, intslice);
    boost::sort::block_indirect_sort(intslice.begin(), intslice.end(), NUM_THREADS);
    auto realnum = std::distance(intslice.begin(), std::lower_bound(intslice.begin(), intslice.end(), DUMMY));
    // move pivots to the end of the memory
    VectorSlice sampleforquantile(intslice, 0, realnum);
    return GetQuantile(sampleforquantile, p0);
}

void OneLevel::Sort(std::vector<int64_t> &input, std::vector<int64_t> &output, SortType sorttype, bool printInfo)
{
    if (printInfo)
    {
        std::string typestr = sorttype == TIGHT ? "Tight" : "Loose";
        std::cout << "------------OneLevel " << typestr << " Sort ------------" << std::endl;
        Tick("Total");
        Tick("Decide pivots");
    }
    auto pivots = GetPivots(input, sorttype);
    if (printInfo)
    {
        Tick("Decide pivots");
        Tick("Partition");
    }
    auto buckets = Partition(input, pivots, true);
    if (printInfo)
    {
        Tick("Partition");
        Tick("Final sorting");
    }
    FinalSorting(buckets, output, sorttype);
    if (printInfo)
    {
        Tick("Final sorting");
        Tick("Total");
        std::cout << "Num IOs: " << (float)m_iom.GetNumIOs() * B / N << "*N/B" << std::endl;
        m_iom.ClearIO();
    }
    for (auto vs : buckets)
        delete vs;
}
