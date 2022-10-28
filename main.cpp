#include <iostream>
#include <vector>
#include "global.h"
#include "memory.h"
#include "ODS.h"
#include <random>
#include <algorithm>
#include <cassert>
#include <chrono>
using namespace std::chrono;

using namespace std;

void print(std::vector<uint64_t> &v, int limit = 100)
{
    if (limit > v.size())
        limit = v.size();
    for (int i = 0; i < limit; i++)
        cout << v[i] << ",";
    cout << endl;
}

void CheckOutput(vector<uint64_t> &output, OneLevel::SortType sorttype)
{
    uint64_t v = 1;
    for (auto t: output)
    {
        if (t == DUMMY)
            assert(sorttype==OneLevel::LOOSE);
        else
            assert(t == v++);
    }
}

// arg1: external number of elements (divided by the max number of elements in the internal memory number), deafult: 128
// arg2: internal memory size (in terms of MB), default: 128
// arg3: block size (in terms of elemens), default: 128

void OneLevelExp(int argc, char *argv[])
{
    uint64_t M, N, B, c, m;
    // N=128*M by default
    c = argc > 1 ? atoll(argv[1]) : 16;
    // M=128MB by default
    m = argc > 2 ? atoll(argv[2]) : 16;
    M = (m<<20)/ sizeof(uint64_t);
    // B=4 by default
    B = argc > 3 ? atoll(argv[3]) : 4;
    N = M * c;
    int sigma = 40;
    Tick("Total");
    IOManager iom (M,B);
    OneLevel ods(iom,N,B,sigma);
    vector<uint64_t> input(N);
    vector<uint64_t> output;
    for(uint64_t i=0;i<N;i++)
        input[i] = N-i;
    //random_shuffle(input.begin(), input.end());
    ods.Sort(input, output, OneLevel::LOOSE);
    CheckOutput(output, OneLevel::LOOSE);
    Tick("Total");
    cout << "Num IOs: " << (float)iom.GetNumIOs()*B/N << "*N/B" << endl;
    iom.ClearIO();
}


int main(int argc, char *argv[])
{
    cout << "Number of threads: " << NUM_THREADS << endl;
    OneLevelExp(argc, argv);
    return 0;
}