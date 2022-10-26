#include <iostream>
#include <vector>
#include "global.h"
#include "memory.h"
#include "ODS.h"
#include "test.h"

using namespace std;


void ReadArgs(int argc, char *argv[])
{
    uint64_t M, N, B, c;
    // N=128*M by default
    c = argc > 1 ? atoll(argv[1]): 128;
    // M=128MB by default
    M = argc > 2 ? atoll(argv[2]): (128<<20/sizeof(uint64_t));
    // B=4 by default
    B = argc > 3 ? atoll(argv[3]): 4;
    N = M * c;
}
int main(int argc, char *argv[])
{
    TestQuantile();
    return 0;
}