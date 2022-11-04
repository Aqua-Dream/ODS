#pragma once
#include <cstdint>
#include <vector>
#include <stdexcept>

template <typename T>
class VectorSlice
{
public:
    // fixRange: if end index out of range, fix it
    VectorSlice(std::vector<T> & source, int64_t ifirst, int64_t length, bool fixRange=false);
    VectorSlice(std::vector<T> & source, typename std::vector<T>::iterator itbegin, typename std::vector<T>::iterator itend);
    VectorSlice(VectorSlice & source, int64_t ifirst, int64_t length, bool fixRange=false);
    VectorSlice(VectorSlice & source, typename std::vector<T>::iterator itbegin, typename std::vector<T>::iterator itend);
    int64_t size();
    int64_t GetStartPos();
    int64_t GetEndPos();
    T &operator[]( int64_t index );
    // fillDummy: if size not full used, fill with dummy
    void CopyDataFrom(VectorSlice &vs, bool fillDummy=false);
    typename std::vector<T>::iterator begin();
    typename std::vector<T>::iterator end();
private:
    int64_t m_ifirst;
    int64_t m_length;
    std::vector<T>& m_source;
};

template <typename T>
class IOManager
{
public:
    IOManager(int64_t internalMemorySize, int64_t blockSize);
    void DataTransfer(VectorSlice<T> &from, VectorSlice<T> &to);
    int64_t GetNumIOs();
    std::vector<T>& GetInternalMemory();
    void PrintInternalMemory(int limit=100);
    void ClearIO();
    int64_t m_numIOs;
private:
    int64_t m_blocksize;
    std::vector<T> m_intmem;
};