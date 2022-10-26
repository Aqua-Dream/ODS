#pragma once
#include <cstdint>
#include <vector>
#include <stdexcept>

class VectorSlice
{
public:
    // fixRange: if end index out of range, fix it
    VectorSlice(std::vector<uint64_t> & source, uint64_t ifirst, uint64_t length, bool fixRange=false);
    VectorSlice(VectorSlice & source, uint64_t ifirst, uint64_t length, bool fixRange=false);
    uint64_t size();
    uint64_t GetStartPos();
    uint64_t GetEndPos();
    uint64_t &operator[]( uint64_t index );
    // fillDummy: if size not full used, fill with dummy
    void CopyDataFrom(VectorSlice &vs, bool fillDummy=false);
    std::vector<uint64_t>::iterator begin();
    std::vector<uint64_t>::iterator end();
private:
    uint64_t m_ifirst;
    uint64_t m_length;
    std::vector<uint64_t>& m_source;
};

class IOManager
{
public:
    IOManager(uint64_t internalMemorySize, uint64_t blockSize);
    void DataTransfer(VectorSlice &from, VectorSlice &to);
    uint64_t GetNumIOs();
    std::vector<uint64_t>& GetInternalMemory();
    void print(int limit=100);
private:
    uint64_t m_numIOs;
    uint64_t m_blocksize;
    std::vector<uint64_t> m_intmem;
};
