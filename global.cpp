#include <random>
#include <cstdint>
#include <vector>
#include "global.h"
#include <chrono>
#include <iostream>
#include <unordered_map>

int NUM_THREADS;
std::random_device dev;
std::mt19937 rng(dev());
std::unordered_map<std::string, std::chrono::system_clock::time_point> tick_table;

// "end" not included
int64_t RandRange(int64_t start, int64_t end)
{
    std::uniform_int_distribution<int64_t> distr(start, end - 1);
    return distr(rng);
}

void Tick(std::string name)
{
    using namespace std::chrono;
    auto it = tick_table.find(name);
    if (it != tick_table.end())
    {
        auto duration = system_clock::now() - it->second;
        tick_table.erase(it);
        int64_t elaspe = duration_cast<milliseconds>(duration).count();
        std::cout << name << ": " << elaspe << " ms" << std::endl;
    }
    else
        tick_table[name] = std::chrono::system_clock::now();
}

int64_t ceil_divide(int64_t n, int64_t q)
{
    int64_t p = n / q;
    if (n % q == 0)
        return p;
    return p + 1;
}