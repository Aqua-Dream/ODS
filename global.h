#pragma once
#include <cstdint>
#include <string>

// "end" not included
uint64_t RandRange(uint64_t start, uint64_t end);
void Tick(std::string name);

const uint64_t DUMMY = 0xCCCCCCCC; 