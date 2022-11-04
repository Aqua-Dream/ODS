#include "memory.h"
#include "global.h"
#include <cmath>
#include <iostream>

template <typename T>
VectorSlice<T>::VectorSlice(std::vector<T> &source, int64_t ifirst, int64_t length, bool fixRange)
    : m_ifirst{ifirst}, m_length{length}, m_source{source}
{
    if (source.size() < ifirst + length)
    {
        if (fixRange)
            this->m_length = source.size() - ifirst;
        else
            throw std::invalid_argument("Index out of range when creating a slice!");
    }
}

template <typename T>
VectorSlice<T>::VectorSlice(std::vector<T> &source, typename std::vector<T>::iterator itbegin, typename std::vector<T>::iterator itend)
    : m_source{source}
{
    m_ifirst = std::distance(m_source.begin(), itbegin);
    m_length = std::distance(itbegin, itend);
}

template <typename T>
VectorSlice<T>::VectorSlice(VectorSlice &source, int64_t ifirst, int64_t length, bool fixRange)
    : m_ifirst{ifirst + source.m_ifirst}, m_length{length}, m_source{source.m_source}
{
    if (source.size() < ifirst + length)
    {
        if (fixRange)
            this->m_length = source.size() - ifirst;
        else
            throw std::invalid_argument("Index out of range when creating a slice!");
    }
}

template <typename T>
VectorSlice<T>::VectorSlice(VectorSlice &source, typename std::vector<T>::iterator itbegin, typename std::vector<T>::iterator itend)
    : m_source{source.m_source}
{
    m_ifirst = std::distance(m_source.begin(), itbegin);
    m_length = std::distance(itbegin, itend);
}

template <typename T>
int64_t VectorSlice<T>::size() { return m_length; }

template <typename T>
int64_t VectorSlice<T>::GetStartPos() { return m_ifirst; }

template <typename T>
int64_t VectorSlice<T>::GetEndPos() { return m_ifirst + m_length; }

template <typename T>
T &VectorSlice<T>::operator[](int64_t index)
{
    if (this->m_length <= index)
        throw std::invalid_argument("Index out of range when accessing a slice!");
    return m_source[m_ifirst + index];
}

template <typename T>
typename std::vector<T>::iterator VectorSlice<T>::begin() { return m_source.begin() + m_ifirst; }

template <typename T>
typename std::vector<T>::iterator VectorSlice<T>::end() { return m_source.begin() + m_ifirst + m_length; }

// fixRange: if end index out of range, fix
template <typename T>
void VectorSlice<T>::CopyDataFrom(VectorSlice &vs, bool fillDummy)
{
    if ((!fillDummy && vs.m_length != this->m_length) || vs.m_length > this->m_length)
        throw std::invalid_argument("Index out of range when copying a slice!");
    auto copyEndIt = std::copy(vs.begin(), vs.end(), this->begin());
    std::fill(copyEndIt, this->end(), DUMMY<T>());
}

template <typename T>
IOManager<T>::IOManager(int64_t memory_size, int64_t block_size)
{
    m_blocksize = block_size;
    m_intmem.resize(memory_size);
    m_numIOs = 0;
}

template <typename T>
void IOManager<T>::DataTransfer(VectorSlice<T> &from, VectorSlice<T> &to)
{
    to.CopyDataFrom(from, true);
    m_numIOs += (int)ceil(to.size() / m_blocksize);
}

template <typename T>
int64_t IOManager<T>::GetNumIOs() { return m_numIOs; }

template <typename T>
void IOManager<T>::PrintInternalMemory(int limit)
{
    if (limit > m_intmem.size())
        limit = m_intmem.size();
    for (int i = 0; i < limit - 1; i++)
        std::cout << m_intmem[i] << ", ";
    std::cout << m_intmem[limit - 1] << std::endl;
}

template <typename T>
void IOManager<T>::ClearIO()
{
    m_numIOs = 0;
}

template <typename T>
std::vector<T> &IOManager<T>::GetInternalMemory() { return m_intmem; }

template class VectorSlice<int32_t>;
template class IOManager<int32_t>;

template class VectorSlice<int64_t>;
template class IOManager<int64_t>;