#include "ODS.h"
#include <cmath>
#include <stdexcept>
#include <boost/sort/sort.hpp>
#include <string>

TwoLevel::TwoLevel(IOManager &iom, int64_t dataSize, int blockSize, int sigma)
    : ObliDistSort(iom, dataSize, blockSize, sigma)
{
    num_levels = 2;
    double kappa = (sigma + 1) * 0.693147;
    // has a solution iff a<0.012858
    double a = 2.0 * B * sqrt(N / M) * (kappa + 1 + 1.5 * log(N / M)) / M;
    beta = argsolver(a);
    if (beta <= 0)
        throw std::invalid_argument("Invalid parameters for two-level sorting!");
    alpha = (kappa + 1 + log(N)) * 4 * (1 + 1 / beta) * (1 / beta + 2) / M;
    uint64_t sample_size = (int64_t)(1.5 * alpha * ceil_divide(N, M) * M);
    while (alpha >= 2 * beta)
    {
        beta *= 2;
        alpha = (kappa + 1 + log(N)) * 4 * (1 + 1 / beta) * (1 / beta + 2) / M;
    }
    if (beta >= 1)
        throw std::invalid_argument("Invalid parameters for two-level sorting!");
    p = p0 = (int)ceil(sqrt((1 + 2 * beta) * N / M));
    if (p0 < 2)
        p0 = 2;
    if (p < 2)
        p = 2;
    int64_t memload = ceil(M / (1 + 2 * beta) / B) * B; // be the multiple of B
    int64_t num_memloads = ceil_divide(N, memload);
    int64_t unit = ceil(float(M) / p0);
    int64_t bucketsize = unit * num_memloads;
    memload = M;
    num_memloads = ceil_divide(bucketsize, memload);
    unit = ceil(float(M) / p);
    bucketsize = unit * num_memloads;
    if (bucketsize > M)
        p++;
}

std::vector<int64_t> TwoLevel::GetPivots(std::vector<int64_t> &data, SortType sorttype)
{
    std::vector<int64_t> sample, sample_sorted;
    Sample(data, sample, sorttype);
    OneLevel subsorter(m_iom, sample.size(), B, sigma + 1);
    subsorter.Sort(sample, sample_sorted, TIGHT);
    VectorSlice vs(sample_sorted, 0, sample_sorted.size());
    int num_total_pivots = p0 * p;
    m_iom.m_numIOs += num_total_pivots;
    return GetQuantile(vs, num_total_pivots);
}

void TwoLevel::Sort(std::vector<int64_t> &input, std::vector<int64_t> &output, SortType sorttype, bool printInfo)
{
    if (printInfo)
    {
        std::string typestr = sorttype == TIGHT ? "Tight" : "Loose";
        std::cout << "------------ TwoLevel " << typestr << " Sort ------------" << std::endl;
        Tick("Total");
        Tick("Decide pivots");
    }
    auto pivots = GetPivots(input, sorttype);
    if (printInfo)
    {
        Tick("Decide pivots");
        Tick("First level partition");
    }
    std::vector<std::vector<int64_t> *> bucketsFinal(p0 * p);
    std::vector<int64_t> pivotsLevelOne(p0 - 1);
    for (int i = 1; i < p0; i++)
        pivotsLevelOne[i - 1] = pivots[i * p - 1];
    auto bucketsLevelOne = Partition(input, pivotsLevelOne, true);
    if (printInfo)
    {
        Tick("First level partition");
        Tick("Second level partition");
    }
    for (int i = 0; i < p0; i++)
    {
        auto bucket = bucketsLevelOne[i];
        std::vector<int64_t> subpivots(pivots.begin() + i * p, pivots.begin() + (i + 1) * p - 1);
        auto subbuckets = Partition(*bucket, subpivots, false);
        std::copy_n(subbuckets.begin(), p, bucketsFinal.begin() + i * p);
    }
    if (printInfo)
    {
        Tick("Second level partition");
        Tick("Final sorting");
    }
    FinalSorting(bucketsFinal, output, sorttype);
    if (printInfo)
    {
        Tick("Final sorting");
        Tick("Total");
        std::cout << "Num IOs: " << (float)m_iom.GetNumIOs() * B / N << "*N/B" << std::endl;
        m_iom.ClearIO();
    }
    for (auto vs : bucketsFinal)
        delete vs;
}