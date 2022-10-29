#include <iostream>
#include <vector>
#include "global.h"
#include "memory.h"
#include "ODS.h"
#include <random>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <boost/program_options.hpp>

using namespace std::chrono;

using namespace std;
namespace po = boost::program_options;

po::variables_map read_options(int argc, char *argv[])
{
    int m, c;
    po::variables_map vm;
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()("help,h", "Show help message")("m,m", po::value<int>()->default_value(16), "Internal memory size (MB)")("c,c", po::value<int>()->default_value(8), "The value N/M")("block_size,B", po::value<uint64_t>()->default_value(4), "Block size (in terms of elements)")("num_threads,T", po::value<int>(&NUM_THREADS)->default_value(4), "Number of threads")("sigma,s", po::value<int>()->default_value(40), "Failure probability upper bound: 2^(-sigma)");
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        if (vm.count("help"))
        {
            cout << desc << endl;
            exit(0);
        }
    }
    catch (exception &e)
    {
        cerr << "error: " << e.what() << endl;
        exit (1);
    }
    catch (...)
    {
        cerr << "Exception of unknown type!\n";
        exit(-1);
    }
    return vm;
}

void CheckOutput(vector<uint64_t> &output, OneLevel::SortType sorttype)
{
    uint64_t v = 1;
    for (auto t : output)
    {
        if (t == DUMMY)
            assert(sorttype == OneLevel::LOOSE);
        else
            assert(t == v++);
    }
}

// arg1: external number of elements (divided by the max number of elements in the internal memory number), deafult: 128
// arg2: internal memory size (in terms of MB), default: 128
// arg3: block size (in terms of elemens), default: 128

void OneLevelExp(po::variables_map vm)
{
    uint64_t M = (vm["m"].as<int>()<<20)/sizeof(uint64_t);
    uint64_t N = vm["c"].as<int>() * M;
    uint64_t B = vm["block_size"].as<uint64_t> ();
    int sigma = vm["sigma"].as<int>();
    Tick("Total");
    IOManager iom(M, B);
    OneLevel ods(iom, N, B, sigma);
    vector<uint64_t> input(N);
    vector<uint64_t> output;
    for (uint64_t i = 0; i < N; i++)
        input[i] = N - i;
    // random_shuffle(input.begin(), input.end());
    ods.Sort(input, output, OneLevel::LOOSE);
    CheckOutput(output, OneLevel::LOOSE);
    Tick("Total");
    cout << "Num IOs: " << (float)iom.GetNumIOs() * B / N << "*N/B" << endl;
    iom.ClearIO();
}

int main(int argc, char *argv[])
{
    auto vm = read_options(argc, argv);
    cout << "Number of threads: " << NUM_THREADS << endl;
    OneLevelExp(vm);
    return 0;
}