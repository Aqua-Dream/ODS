#include "ODS.h"
#include <cmath>
#include <stdexcept>
#include <boost/sort/sort.hpp>

template <typename T>
OneLevel<T>::OneLevel(IOManager<T> &iom, int64_t dataSize, int blockSize, int sigma)
    : ObliDistSort<T>(iom, dataSize, blockSize, sigma)
{
    int64_t M = this->M;
    int64_t N  = this->N;
    int B = this->B;
    this->num_levels = 1;
    double kappa = sigma * 0.693147;
    // has a solution iff a<0.012858
    double a = 2.0 * N * B / M * (kappa + 1 + 2 * log(N / M)) / M;
    double beta = argsolver(a);
    if (beta <= 0)
        throw std::invalid_argument("Invalid parameters for one-level sorting!");
    double alpha = (kappa + 1 + log(N)) * 4 * (1 + 1 / beta) * (1 / beta + 2) / M;
    while (this->GetSampleSize() > M || alpha >= 2*beta)
    {
        beta *= 2;
        alpha = (kappa + 1 + log(N)) * 4 * (1 + 1 / beta) * (1 / beta + 2) / M;
    }
    if (beta >= 1)
        throw std::invalid_argument("Invalid parameters for one-level sorting!");
    int p0 = (int)ceil((1 + 2 * beta) * N / M);
    int64_t memload = ceil(M / (1 + 2 * beta) / B) * B; // be the multiple of B
    int64_t num_memloads = ceil_divide(N, memload);
    int64_t unit = ceil(float(M) / p0);
    if (unit * num_memloads > M)
        p0++;
    this->p0 = p0;
    this->alpha = alpha;
    this->beta = beta;
}

template <typename T>
std::vector<T> OneLevel<T>::GetPivots(std::vector<T> &data, SortType sorttype)
{
    std::vector<T> sample;
    this->Sample(data, sample, sorttype);
    if (sample.size() >this-> M)
        throw std::invalid_argument("Sample does not fit into memory!");
    VectorSlice<T> sampleslice(sample, 0, sample.size());
    VectorSlice<T> intslice(this->m_intmem, 0, sample.size());
    this->m_iom.DataTransfer(sampleslice, intslice);
    boost::sort::block_indirect_sort(intslice.begin(), intslice.end(), NUM_THREADS);
    auto realnum = std::distance(intslice.begin(), std::lower_bound(intslice.begin(), intslice.end(), DUMMY<T>()));
    // move pivots to the end of the memory
    VectorSlice<T> sampleforquantile(intslice, 0, realnum);
    return GetQuantile<T>(sampleforquantile, this->p0);
}

template <typename T>
void OneLevel<T>::Sort(std::vector<T> &input, std::vector<T> &output, SortType sorttype, bool printInfo)
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
    auto buckets = this->Partition(input, pivots, true);
    if (printInfo)
    {
        Tick("Partition");
        Tick("Final sorting");
    }
    this->FinalSorting(buckets, output, sorttype);
    if (printInfo)
    {
        Tick("Final sorting");
        Tick("Total");
        std::cout << "Num IOs: " << (float)this->m_iom.GetNumIOs() * this->B / this->N << "*N/B" << std::endl;
        this->m_iom.ClearIO();
    }
    for (auto vs : buckets)
        delete vs;
}

template class OneLevel<int32_t>;
template class OneLevel<int64_t>;