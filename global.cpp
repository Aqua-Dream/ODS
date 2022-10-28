#include <random>
#include <cstdint>
#include <vector>
#include "global.h"
#include <chrono>
#include <iostream>
#include <unordered_map>

std::random_device dev;
std::mt19937 rng(0); 
std::unordered_map<std::string, std::chrono::system_clock::time_point> tick_table;

// "end" not included
uint64_t RandRange(uint64_t start, uint64_t end)
{
    std::uniform_int_distribution<uint64_t>  distr(start, end-1);
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
    tick_table[name] = std::chrono::system_clock::now();
}
