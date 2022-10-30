#pragma once
#include <cstdint>
#include <vector>
#include <stdexcept>

class VectorSlice
{
public:
    // fixRange: if end index out of range, fix it
    VectorSlice(std::vector<int64_t> & source, int64_t ifirst, int64_t length, bool fixRange=false);
    VectorSlice(std::vector<int64_t> & source, std::vector<int64_t>::iterator itbegin, std::vector<int64_t>::iterator itend);
    VectorSlice(VectorSlice & source, int64_t ifirst, int64_t length, bool fixRange=false);
    VectorSlice(VectorSlice & source, std::vector<int64_t>::iterator itbegin, std::vector<int64_t>::iterator itend);
    int64_t size();
    int64_t GetStartPos();
    int64_t GetEndPos();
    int64_t &operator[]( int64_t index );
    // fillDummy: if size not full used, fill with dummy
    void CopyDataFrom(VectorSlice &vs, bool fillDummy=false);
    std::vector<int64_t>::iterator begin();
    std::vector<int64_t>::iterator end();
private:
    int64_t m_ifirst;
    int64_t m_length;
    std::vector<int64_t>& m_source;
};

class IOManager
{
public:
    IOManager(int64_t internalMemorySize, int64_t blockSize);
    void DataTransfer(VectorSlice &from, VectorSlice &to);
    int64_t GetNumIOs();
    std::vector<int64_t>& GetInternalMemory();
    void PrintInternalMemory(int limit=100);
    void ClearIO();
    int64_t m_numIOs;
private:
    int64_t m_blocksize;
    std::vector<int64_t> m_intmem;
};
