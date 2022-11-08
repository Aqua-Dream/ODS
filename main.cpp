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
        desc.add_options()("help,h", "Show help message")("m,m", po::value<int>()->default_value(8), "Internal memory size (MB)")("c,c", po::value<int>()->default_value(16), "The value N/M")("block_size,B", po::value<int>()->default_value(4), "Block size (in terms of elements)")("num_threads,T", po::value<int>(&NUM_THREADS)->default_value(4), "Number of threads")("sigma,s", po::value<int>()->default_value(40), "Failure probability upper bound: 2^(-sigma)")("onelevel", "OneLevel test")("twolevel", "OneLevel test")("failure", "Failure test")("reps,r", po::value<int>()->default_value(100), "Number of repetitions in failure test.")("wide,w", "Use wide data type int64_t instead of int32_t.");
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
        cerr << "Error: " << e.what() << endl;
        exit(1);
    }
    catch (...)
    {
        cerr << "Exception of unknown type!\n";
        exit(-1);
    }
    return vm;
}

template <typename T>
void CheckOutput(vector<T> &output, SortType sorttype)
{
    T v = 1;
    for (auto t : output)
    {
        if (t == DUMMY<T>() && sorttype != LOOSE)
            throw "Dummy found in tight sort!";
        else if (t != DUMMY<T>() && t != v++)
            throw "Value error!";
    }
}

template <typename T>
void OneLevelExp(po::variables_map vm)
{
    int64_t M = ((int64_t)(vm["m"].as<int>()) << 20) / sizeof(T);
    int64_t N = vm["c"].as<int>() * M;
    int B = vm["block_size"].as<int>();
    int sigma = vm["sigma"].as<int>();

    IOManager<T> iom(M, B);
    OneLevel<T> ods(iom, N, B, sigma);
    vector<T> input(N);
    vector<T> output;
#pragma omp parallel for
    for (T i = 0; i < N; i++)
        input[i] = N - i;

    ods.Sort(input, output, TIGHT, true);
    CheckOutput(output, TIGHT);

    ods.Sort(input, output, LOOSE, true);
    CheckOutput(output, LOOSE);
}

template <typename T>
void TwoLevelExp(po::variables_map vm)
{
    int64_t M = ((int64_t)(vm["m"].as<int>()) << 20) / sizeof(T);
    int64_t N = vm["c"].as<int>() * M;
    int B = vm["block_size"].as<int>();
    int sigma = vm["sigma"].as<int>();

    IOManager<T> iom(M, B);
    TwoLevel<T> ods(iom, N, B, sigma);
    vector<T> input(N);
    vector<T> output;
#pragma omp parallel for
    for (T i = 0; i < N; i++)
        input[i] = N - i;

    ods.Sort(input, output, TIGHT, true);
    CheckOutput(output, TIGHT);

    ods.Sort(input, output, LOOSE, true);
    CheckOutput(output, LOOSE);
}

template <typename T>
void FailureTest(po::variables_map vm)
{
    int reps = vm["reps"].as<int>();
    int num_total_fails = 0;
    int64_t M = ((int64_t)(vm["m"].as<int>()) << 20) / sizeof(T);
    int64_t N = vm["c"].as<int>() * M;
    int B = vm["block_size"].as<int>();
    int sigma = vm["sigma"].as<int>();
    const int step = 100;
    vector<T> input(N);
#pragma omp parallel for
    for (int64_t i = 0; i < N; i++)
        input[i] = N - i;
    IOManager<T> iom(M, B);
    OneLevel<T> ods(iom, N, B, sigma);
    vector<T> output;
    reps = ceil_divide(reps, step) * step;
    for (int i = 0; i < reps/step; i++)
    {
        int num_step_fails = 0;
        for (int j = 0; j < step; j++)
        {
            try
            {
                ods.Sort(input, output, LOOSE);
            }
            catch (...)
            {
                num_step_fails++;
            }
        }
        cout << num_step_fails << " fails in the " << i+1 << "-th " << step << " trials." << endl;
        num_total_fails += num_step_fails;
    }
    cout << num_total_fails << " fails out of " << reps << " trials." << endl;
}

using FP = void (*)(po::variables_map);

template <typename T>
vector<FP> GetTestFunctions(po::variables_map vm)
{
    vector<FP> tester;
    if (vm.count("onelevel"))
        tester.push_back(OneLevelExp<T>);
    if (vm.count("twolevel"))
        tester.push_back(TwoLevelExp<T>);
    if (vm.count("failure"))
        tester.push_back(FailureTest<T>);
    return tester;
}

int main(int argc, char *argv[])
{
    auto vm = read_options(argc, argv);
    vector<FP> tester;
    if (vm.count("wide"))
        tester = GetTestFunctions<int64_t>(vm);
    else
        tester = GetTestFunctions<int32_t>(vm);
    for (auto func : tester)
    {
        try
        {
            func(vm);
        }
        catch (exception &e)
        {
            cerr << "Error: " << e.what() << endl;
        }
        catch (...)
        {
            cerr << "Unknown error occurs." << endl;
        }
    }
    return 0;
}