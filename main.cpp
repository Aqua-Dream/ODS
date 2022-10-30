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
#include <omp.h>

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
        desc.add_options()("help,h", "Show help message")("m,m", po::value<int>()->default_value(8), "Internal memory size (MB)")("c,c", po::value<int>()->default_value(16), "The value N/M")("block_size,B", po::value<int>()->default_value(4), "Block size (in terms of elements)")("num_threads,T", po::value<int>(&NUM_THREADS)->default_value(4), "Number of threads")("sigma,s", po::value<int>()->default_value(40), "Failure probability upper bound: 2^(-sigma)")("onelevel", "OneLevel test")("twolevel", "OneLevel test")("failure", "Failure test")("reps,r", po::value<int>()->default_value(100), "Number of repetitions in failure test.");
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        if (vm.count("help"))
        {
            cout << desc << endl;
            exit(0);
        }
        omp_set_num_threads(NUM_THREADS);
    }
    catch (exception &e)
    {
        cerr << "error: " << e.what() << endl;
        exit(1);
    }
    catch (...)
    {
        cerr << "Exception of unknown type!\n";
        exit(-1);
    }
    return vm;
}

void CheckOutput(vector<int64_t> &output, OneLevel::SortType sorttype)
{
    int64_t v = 1;
    for (auto t : output)
    {
        if (t == DUMMY && sorttype != OneLevel::LOOSE)
            throw "Dummy found in tight sort!";
        else if (t != DUMMY && t != v++)
            throw "Value error!";
    }
}

void OneLevelExp(po::variables_map vm)
{
    int64_t M = (vm["m"].as<int>() << 20) / sizeof(int64_t);
    int64_t N = vm["c"].as<int>() * M;
    int B = vm["block_size"].as<int>();
    int sigma = vm["sigma"].as<int>();

    IOManager iom(M, B);
    OneLevel ods(iom, N, B, sigma);
    vector<int64_t> input(N);
    vector<int64_t> output;
#pragma omp parallel for
    for (int64_t i = 0; i < N; i++)
        input[i] = N - i;

    ods.Sort(input, output, OneLevel::TIGHT, true);
    CheckOutput(output, OneLevel::TIGHT);

    ods.Sort(input, output, OneLevel::LOOSE, true);
    CheckOutput(output, OneLevel::LOOSE);
}

void TwoLevelExp(po::variables_map vm)
{
    int64_t M = (vm["m"].as<int>() << 20) / sizeof(int64_t);
    int64_t N = vm["c"].as<int>() * M;
    int B = vm["block_size"].as<int>();
    int sigma = vm["sigma"].as<int>();

    IOManager iom(M, B);
    TwoLevel ods(iom, N, B, sigma);
    vector<int64_t> input(N);
    vector<int64_t> output;
#pragma omp parallel for
    for (int64_t i = 0; i < N; i++)
        input[i] = N - i;

    ods.Sort(input, output, TwoLevel::TIGHT, true);
    CheckOutput(output, TwoLevel::TIGHT);

    ods.Sort(input, output, TwoLevel::LOOSE, true);
    CheckOutput(output, TwoLevel::LOOSE);
}

void FailureTest(po::variables_map vm)
{
    int reps = vm["reps"].as<int>();
    int num_fails = 0;
#pragma omp parallel for reduction(+:num_fails)
    for (int i = 0; i < reps; i++)
    {
        try
        {
            int64_t M = (vm["m"].as<int>() << 20) / sizeof(int64_t);
            int64_t N = vm["c"].as<int>() * M;
            int B = vm["block_size"].as<int>();
            int sigma = vm["sigma"].as<int>();
            IOManager iom(M, B);
            OneLevel ods(iom, N, B, sigma);
            vector<int64_t> input(N);
            vector<int64_t> output;
            for (int64_t i = 0; i < N; i++)
                input[i] = N - i;
            ods.Sort(input, output, OneLevel::LOOSE);
        }
        catch(...)
        {
            num_fails++;
        }
    }
    cout << num_fails << " fails out of " << reps << " trials." << endl;
}

int main(int argc, char *argv[])
{
    auto vm = read_options(argc, argv);
    if (vm.count("onelevel"))
        OneLevelExp(vm);
    if (vm.count("twolevel"))
        TwoLevelExp(vm);
    if (vm.count("failure"))
        FailureTest(vm);
    return 0;
}