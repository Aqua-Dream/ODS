#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "Feistel.h"
#include "global.h"
#include "memory.h"
#include "sampling.h"
#include "ODS.h"
#include "argsolver.h"

void TestSampling()
{
    int M=10,B=2,N=50,n=8;
    IOManager iom(M, B);  
    vector<uint64_t> extmem(N);
    auto& intmem = iom.GetInternalMemory();
    for(int i=0;i<N;i++)
    {
        extmem[i] = RandRange(1,100);
        cout << extmem[i] << ",";
    }
    cout << endl;
    SmallSample(iom, extmem, B, n, false);
    for(int i=0;i<M;i++)
    {
        cout << intmem[i] << ",";
    }
    cout << endl;
    cout << iom.GetNumIOs() << endl;
}

void TestInternalPartition()
{
    int M=30,B=2,N=100,n=16;
    IOManager iom(M, B);
    auto& intmem = iom.GetInternalMemory();
    for(int i=0;i<20;i++)
    {
        intmem[i] = RandRange(1,100);
        if(intmem[i]%5==0)
            intmem[i]++;
        cout << intmem[i] << ",";
    }
    cout << endl;
    for(int i=20;i<25;i++)
    {
        intmem[i] = (i-20)*5+90;
        cout << intmem[i] << ",";
    }
    cout << endl;
    for(int i=25;i<30;i++)
    {
        //intmem[i] = (i-25)*20;
        cout << intmem[i] << ",";
    }
    cout << endl;
    VectorSlice data(intmem, 0, 20);
    VectorSlice pivots(intmem, 20, 5);
    VectorSlice posList(intmem, 25, 5);
    InternalPartition(data, pivots, posList);
    for(int i=0;i<20;i++)
    {
        //intmem[i] = RandRange(1,100);
        cout << intmem[i] << ",";
    }
    cout << endl;
    for(int i=20;i<25;i++)
    {
        //intmem[i] = (i-20)*15;
        cout << intmem[i] << ",";
    }
    cout << endl;
    for(int i=25;i<30;i++)
    {
        //intmem[i] = (i-25)*20;
        cout << intmem[i] << ",";
    }
    cout << endl;

}



// ans by python: 0.011699, 0.0172392, 0.025643, 0.0387563, 0.060241
void TestArgSolver()
{
    uint64_t M = 16777216;
    uint64_t B = 4;
    double kappa = 27.8;
    for(uint64_t c=8; c<=512; c*=2)
    {
        double a = 2.0*sqrt((double)c)*B*(kappa+2+1.5*log(c))/M;
        double ans = argsolver(a);
        cout << ans << ",";
    }
}


void TestQuantile()
{
    int N=50, q=4;
    vector<uint64_t> data(N);
    for(int i=0;i<N;i++)
    {
        data[i] = i+1;
    }
    Quantile(data, N, q);
    for (auto v : data)
        cout << v << ",";
    cout << endl;
    for(int i=0;i<N;i++)
    {
        data[i] = i+1;
    }
    Quantile(data, N-1, q+1);
    for (auto v : data)
        cout << v << ",";
    cout << endl;
}
