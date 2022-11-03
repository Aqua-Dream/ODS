#pragma once
#include <cstdint>
#include <string>

// "end" not included
int64_t RandRange(int64_t start, int64_t end);
void Tick(std::string name);

const int64_t DUMMY = INT64_MAX; 
extern int NUM_THREADS;

int64_t ceil_divide(int64_t n, int64_t q);