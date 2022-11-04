#pragma once
#include <cstdint>
#include <string>
#include <limits>

// "end" not included
int64_t RandRange(int64_t start, int64_t end);
void Tick(std::string name);

template<typename T>
constexpr T DUMMY()
{
    return std::numeric_limits<T>::max();
} 

extern int NUM_THREADS;

int64_t ceil_divide(int64_t n, int64_t q);