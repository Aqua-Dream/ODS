#include "memory.h"
#include "global.h"
#include <cmath>
#include <iostream>
#include "omp.h"

VectorSlice::VectorSlice(std::vector<uint64_t> &source, uint64_t ifirst, uint64_t length, bool fixRange)
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

VectorSlice::VectorSlice(VectorSlice &source, uint64_t ifirst, uint64_t length, bool fixRange)
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

uint64_t VectorSlice::size() { return m_length; }
uint64_t VectorSlice::GetStartPos() { return m_ifirst; }
uint64_t VectorSlice::GetEndPos() { return m_ifirst + m_length; }

uint64_t &VectorSlice::operator[](uint64_t index)
{
    if (this->m_length <= index)
        throw std::invalid_argument("Index out of range when accessing a slice!");
    return m_source[m_ifirst + index];
}

std::vector<uint64_t>::iterator VectorSlice::begin() { return m_source.begin() + m_ifirst; }

std::vector<uint64_t>::iterator VectorSlice::end() { return m_source.begin() + m_ifirst + m_length; }

// fixRange: if end index out of range, fix
void VectorSlice::CopyDataFrom(VectorSlice &vs, bool fillDummy)
{
    if ((!fillDummy && vs.m_length != this->m_length) || vs.m_length > this->m_length)
        throw std::invalid_argument("Index out of range when copying a slice!");
    auto copyEndIt = std::copy(vs.begin(), vs.end(), this->begin());
    std::fill(copyEndIt, this->end(), DUMMY);
}

IOManager::IOManager(uint64_t memory_size, uint64_t block_size)
{
    this->m_blocksize = block_size;
    this->m_intmem.resize(memory_size);
    this->m_numIOs = 0;
}

void IOManager::DataTransfer(VectorSlice &from, VectorSlice &to)
{
    to.CopyDataFrom(from, true);
    m_numIOs += (int)ceil(to.size() / m_blocksize);
}

uint64_t IOManager::GetNumIOs() { return m_numIOs; }
void IOManager::PrintInternalMemory(int limit)
{
    if (limit > m_intmem.size())
        limit = m_intmem.size();
    for (int i = 0; i < limit - 1; i++)
        std::cout << m_intmem[i] << ", ";
    std::cout << m_intmem[limit - 1] << std::endl;
}

void IOManager::ClearIO()
{
    m_numIOs = 0;
}

std::vector<uint64_t> &IOManager::GetInternalMemory() { return m_intmem; }