#include "ODS.h"
#include <cmath>
#include <stdexcept>
#include "sampling.h"
#include "argsolver.h"
#include <algorithm>
#include <cassert>

// compute the q-quantile of the first n elements in data, and put them at the end of data
void Quantile(std::vector<uint64_t>& data, uint64_t n, int q)
{
    uint64_t M = data.size();
    uint64_t step = n/q;
    int remain = n%q;
    uint64_t pivotId = M-1;
    // dataIdx<n<=M<int_max, so there is no overflow
    int dataIdx = n;
    int i=0;
    for(;i<remain;i++)
    {
        dataIdx -= step+1; 
        data[pivotId--] = data[dataIdx];
    }
    for(; i<q-1 ;i++)
    {
        dataIdx -= step; 
        data[pivotId--] = data[dataIdx];
    }
    assert(dataIdx == step);
    assert(pivotId == M-q);
}

// return: the number of elements smaller than pivot
uint64_t InternalPartitionLevel(VectorSlice &data, uint64_t pivot)
{
    uint64_t i = 0;
    for (uint64_t j = 0; j < data.size(); ++j)
    {
        if (pivot > data[j])
        {
            uint64_t tmp = data[i];
            data[i] = data[j];
            data[j] = tmp;
            i++;
        }
    }
    return i;
}

void InternalPartition(VectorSlice &data, VectorSlice &pivots, VectorSlice &posList)
{
    uint64_t pivotIdx, pivot, dataIdx;
    pivotIdx = pivots.size() >> 1;
    pivot = pivots[pivotIdx];
    dataIdx = InternalPartitionLevel(data, pivot);
    posList[pivotIdx] = data.GetStartPos() + dataIdx;
    if (dataIdx > 0 && pivotIdx > 0)
    {
        VectorSlice dataLeft(data, 0, dataIdx);
        VectorSlice pivotsLeft(pivots, 0, pivotIdx);
        VectorSlice posLeft(posList, 0, pivotIdx);
        InternalPartition(dataLeft, pivotsLeft, posLeft);
    }
    else if (dataIdx <= 0)
    {
        for (uint64_t i = 0; i < pivotIdx; i++)
            posList[i] = data.GetStartPos();
    }
    if (dataIdx < data.size() && pivotIdx + 1 < pivots.size())
    {
        VectorSlice dataRight(data, dataIdx, data.size() - dataIdx);
        VectorSlice pivotsRight(pivots, pivotIdx + 1, pivots.size() - pivotIdx - 1);
        VectorSlice posRight(posList, pivotIdx + 1, pivots.size() - pivotIdx - 1);
        InternalPartition(dataRight, pivotsRight, posRight);
    }
    else if (dataIdx >= data.size())
    {
        for (uint64_t i = pivotIdx + 1; i < posList.size(); i++)
            posList[i] = data.GetEndPos();
    }
}

OneLevel::OneLevel(IOManager &iom, std::vector<uint64_t> &extdata, uint64_t blockSize, int sigma)
    : N{extdata.size()}, M{iom.GetInternalMemory().size()}, B{blockSize}, m_iom{iom}, m_intmem{iom.GetInternalMemory()}, m_extdata{extdata}
{
    double kappa = sigma * 0.693147;
    double a = 2.0 * N * B / M * (kappa + 1 + 2 * log(N / M)) / M;
    double beta = argsolver(a);
    if (beta == -1)
        throw std::invalid_argument("Invalid parameters for one-level sorting!");
    this->beta = beta;
    this->alpha = (kappa + 1 + log(N)) * 4 * (1 + 1 / beta) * (1 / beta + 2) / M;
    this->p = (int)ceil((1 + 2 * beta) * N / M);
}

void OneLevel::DecidePivots(Type m_type)
{
    uint64_t n = (int)ceil(alpha * N);
    SmallSample(m_iom, m_extdata, B, n, m_type == TIGHT);
    std::sort(m_intmem.begin(), m_intmem.begin() + n);
    // move pivots to the end of the memory
    Quantile(m_intmem, n, p);
}

// std::pair<int, int> OneLevelPartition(int inStructureId, int inSize, std::vector<int> &samples, int sampleSize, int p, int outStructureId1, int is_rec)
// {
//   if (inSize <= M)
//   {
//     return {inSize, 1};
//   }
//   std::cout << "In Onelevel Partition\n";
//   double beta = (!is_rec) ? BETA : _BETA;
//   int hatN = ceil(1.0 * (1 + 2 * beta) * inSize);
//   int M_prime = ceil(1.0 * M / (1 + 2 * beta));
//   int r = ceil(1.0 * log(hatN / M) / log(p));
//   int p0 = ceil(1.0 * hatN / (M * pow(p, r - 1)));
//   quantileCal(samples, 0, sampleSize, p0);
//   int boundary1 = ceil(1.0 * inSize / M_prime);
//   int boundary2 = ceil(1.0 * M_prime / BLOCK_DATA_SIZE);
//   int dataBoundary = boundary2 * BLOCK_DATA_SIZE;
//   int smallSectionSize = M / p0;
//   int bucketSize0 = boundary1 * smallSectionSize;
//   freeAllocate(outStructureId1, outStructureId1, boundary1 * smallSectionSize * p0);

//   int Msize1, Msize2, index1, index2, writeBackNum;
//   int total_blocks = ceil(1.0 * inSize / BLOCK_DATA_SIZE);
//   int *trustedM3 = (int *)malloc(sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
//   memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
//   int *shuffleB = (int *)malloc(sizeof(int) * BLOCK_DATA_SIZE);
//   std::vector<int> partitionIdx;
//   // Finish FFSEM implementation in c++
//   pseudo_init(total_blocks);
//   int index_range = max_num;
//   int k = 0, read_index;
//   for (int i = 0; i < boundary1; ++i)
//   {
//     for (int j = 0; j < boundary2; ++j)
//     {
//       read_index = encrypt(k);
//       while (read_index >= total_blocks)
//       {
//         k += 1;
//         if (k == index_range)
//         {
//           k = -1;
//           break;
//         }
//         read_index = encrypt(k);
//       }
//       if (k == -1)
//       {
//         break;
//       }
//       Msize1 = std::min(BLOCK_DATA_SIZE, inSize - read_index * BLOCK_DATA_SIZE);
//       opOneLinearScanBlock(read_index * BLOCK_DATA_SIZE, &trustedM3[j * BLOCK_DATA_SIZE], Msize1, inStructureId, 0, 0);
//       k += 1;
//       if (k == index_range)
//       {
//         break;
//       }
//     }
//     int blockNum = moveDummy(trustedM3, dataBoundary);
//     quickSortMulti(trustedM3, 0, blockNum - 1, samples, 1, p0, partitionIdx);
//     sort(partitionIdx.begin(), partitionIdx.end());
//     partitionIdx.insert(partitionIdx.begin(), -1);
//     for (int j = 0; j < p0; ++j)
//     {
//       index1 = partitionIdx[j] + 1;
//       index2 = partitionIdx[j + 1];
//       writeBackNum = index2 - index1 + 1;
//       if (writeBackNum > smallSectionSize)
//       {
//         printf("Overflow in small section M/p0: %d", writeBackNum);
//       }
//       opOneLinearScanBlock(j * bucketSize0 + i * smallSectionSize, &trustedM3[index1], writeBackNum, outStructureId1, 1, smallSectionSize - writeBackNum);
//     }
//     memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
//     partitionIdx.clear();
//   }
//   free(trustedM3);
//   free(shuffleB);
//   if (bucketSize0 > M)
//   {
//     printf("Each section size is greater than M, adjst parameters: %d, %d", bucketSize0, M);
//   }
//   return {bucketSize0, p0};
// }

// // TODO: Add TwoLevelPartition
// std::pair<int, int> TwoLevelPartition(int inStructureId, std::vector<std::vector<int>> &pivots, int p, int outStructureId1, int outStructureId2)
// {
//   int M_prime = ceil(1.0 * M / (1 + 2 * BETA));
//   int p0 = p;
//   int boundary1 = ceil(1.0 * N / M_prime);
//   int boundary2 = ceil(1.0 * M_prime / BLOCK_DATA_SIZE);
//   int dataBoundary = boundary2 * BLOCK_DATA_SIZE;
//   int smallSectionSize = M / p0;
//   int bucketSize0 = boundary1 * smallSectionSize;
//   freeAllocate(outStructureId1, outStructureId1, boundary1 * smallSectionSize * p0);
//   int k, Msize1, Msize2, index1, index2, writeBackNum;
//   int blocks_done = 0;
//   int total_blocks = ceil(1.0 * N / BLOCK_DATA_SIZE);
//   int *trustedM3 = (int *)malloc(sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
//   memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
//   int *shuffleB = (int *)malloc(sizeof(int) * BLOCK_DATA_SIZE);
//   std::vector<int> partitionIdx;
//   for (int i = 0; i < boundary1; ++i)
//   {
//     for (int j = 0; j < boundary2; ++j)
//     {
//       if (total_blocks - 1 - blocks_done == 0)
//       {
//         k = 0;
//       }
//       else
//       {
//         k = rand() % (total_blocks - blocks_done);
//       }
//       Msize1 = std::min(BLOCK_DATA_SIZE, N - k * BLOCK_DATA_SIZE);
//       opOneLinearScanBlock(k * BLOCK_DATA_SIZE, &trustedM3[j * BLOCK_DATA_SIZE], Msize1, inStructureId, 0, 0);
//       memset(shuffleB, DUMMY, sizeof(int) * BLOCK_DATA_SIZE);
//       Msize2 = std::min(BLOCK_DATA_SIZE, N - (total_blocks - 1 - blocks_done) * BLOCK_DATA_SIZE);
//       opOneLinearScanBlock((total_blocks - 1 - blocks_done) * BLOCK_DATA_SIZE, shuffleB, Msize2, inStructureId, 0, 0);
//       opOneLinearScanBlock(k * BLOCK_DATA_SIZE, shuffleB, BLOCK_DATA_SIZE, inStructureId, 1, 0);
//       blocks_done += 1;
//       if (blocks_done == total_blocks)
//       {
//         break;
//       }
//     }
//     int blockNum = moveDummy(trustedM3, dataBoundary);
//     quickSortMulti(trustedM3, 0, blockNum - 1, pivots[0], 1, p0, partitionIdx);
//     sort(partitionIdx.begin(), partitionIdx.end());
//     partitionIdx.insert(partitionIdx.begin(), -1);
//     for (int j = 0; j < p0; ++j)
//     {
//       index1 = partitionIdx[j] + 1;
//       index2 = partitionIdx[j + 1];
//       writeBackNum = index2 - index1 + 1;
//       if (writeBackNum > smallSectionSize)
//       {
//         std::cout << "Overflow in small section M/p0: " << writeBackNum << std::endl;
//       }
//       opOneLinearScanBlock(j * bucketSize0 + i * smallSectionSize, &trustedM3[index1], writeBackNum, outStructureId1, 1, smallSectionSize - writeBackNum);
//     }
//     memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
//     partitionIdx.clear();
//   }
//   free(trustedM3);
//   free(shuffleB);
//   // Level2
//   int p1 = p0 * p, readSize, readSize2, k1, k2;
//   int boundary3 = ceil(1.0 * bucketSize0 / M);
//   int bucketSize1 = boundary3 * smallSectionSize;
//   freeAllocate(outStructureId2, outStructureId2, boundary3 * smallSectionSize * p0 * p);
//   // std::vector<int> trustedM1;
//   int *trustedM2 = (int *)malloc(sizeof(int) * M);
//   // TODO: Change memory excceeds? use &trustedM2[M-1-p]
//   int *trustedM2_part = (int *)malloc(sizeof(int) * (p + 1));
//   for (int j = 0; j < p0; ++j)
//   {
//     // trustedM1 = pivots[1+j];
//     for (int k = 0; k < boundary3; ++k)
//     {
//       Msize1 = std::min(M, bucketSize0 - k * M);
//       readSize = (Msize1 < (p + 1)) ? Msize1 : (Msize1 - (p + 1));
//       opOneLinearScanBlock(j * bucketSize0 + k * M, trustedM2, readSize, outStructureId1, 0, 0);
//       k1 = moveDummy(trustedM2, readSize);
//       readSize2 = (Msize1 < (p + 1)) ? 0 : (p + 1);
//       opOneLinearScanBlock(j * bucketSize0 + k * M + readSize, trustedM2_part, readSize2, outStructureId1, 0, 0);
//       k2 = moveDummy(trustedM2_part, readSize2);
//       memcpy(&trustedM2[k1], trustedM2_part, sizeof(int) * k2);
//       quickSortMulti(trustedM2, 0, k1 + k2 - 1, pivots[1 + j], 1, p, partitionIdx);
//       sort(partitionIdx.begin(), partitionIdx.end());
//       partitionIdx.insert(partitionIdx.begin(), -1);
//       for (int ll = 0; ll < p; ++ll)
//       {
//         index1 = partitionIdx[ll] + 1;
//         index2 = partitionIdx[ll + 1];
//         writeBackNum = index2 - index1 + 1;
//         if (writeBackNum > smallSectionSize)
//         {
//           std::cout << "Overflow in small section M/p0: " << writeBackNum << std::endl;
//         }
//         opOneLinearScanBlock((j * p0 + ll) * bucketSize1 + k * smallSectionSize, &trustedM2[index1], writeBackNum, outStructureId2, 1, smallSectionSize - writeBackNum);
//       }
//       memset(trustedM2, DUMMY, sizeof(int) * M);
//       partitionIdx.clear();
//     }
//   }
//   if (bucketSize1 > M)
//   {
//     std::cout << "Each section size is greater than M, adjust parameters: " << bucketSize1 << std::endl;
//   }
//   return {bucketSize1, p1};
// }

// int ObliviousTightSort(int inStructureId, int inSize, int outStructureId1, int outStructureId2)
// {
//   int *trustedM;
//   std::cout << "In ObliviousTightSort\n";
//   if (inSize <= M)
//   {
//     trustedM = (int *)malloc(sizeof(int) * M);
//     opOneLinearScanBlock(0, trustedM, inSize, inStructureId, 0, 0);
//     quickSort(trustedM, 0, inSize - 1);
//     freeAllocate(outStructureId1, outStructureId1, inSize);
//     opOneLinearScanBlock(0, trustedM, inSize, outStructureId1, 1, 0);
//     free(trustedM);
//     return outStructureId1;
//   }
//   std::vector<int> trustedM2;
//   int realNum = Sample(inStructureId, inSize, trustedM2, is_tight);
//   std::pair<int, int> section = OneLevelPartition(inStructureId, inSize, trustedM2, realNum, P, outStructureId1, 0);
//   int sectionSize = section.first;
//   int sectionNum = section.second;
//   // TODO: IN order to reduce memory, can replace outStructureId2 with inStructureId
//   freeAllocate(outStructureId2, outStructureId2, inSize);
//   trustedM = (int *)malloc(sizeof(int) * M);
//   int j = 0, k;
//   std::cout << "In final\n";
//   for (int i = 0; i < sectionNum; ++i)
//   {
//     opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId1, 0, 0);
//     k = moveDummy(trustedM, sectionSize);
//     quickSort(trustedM, 0, k - 1);
//     opOneLinearScanBlock(j, trustedM, k, outStructureId2, 1, 0);
//     j += k;
//     if (j > inSize)
//     {
//       std::cout << "Final error" << std::endl;
//     }
//   }
//   free(trustedM);
//   return outStructureId2;
// }

// // TODO: TightSort2
// int ObliviousTightSort2(int inStructureId, int inSize, int sampleId, int sortedSampleId, int outStructureId1, int outStructureId2)
// {
//   std::cout << "In ObliviousTightSort2 && In SampleRec\n";
//   std::vector<std::vector<int>> pivots;
//   SampleRec(inStructureId, sampleId, sortedSampleId, 1, pivots);
//   std::cout << "In TwoLevelPartition\n";
//   std::pair<int, int> section = TwoLevelPartition(inStructureId, pivots, P, outStructureId1, outStructureId2);
//   std::cout << "Till Partition IOcost: " << IOcost / N * BLOCK_DATA_SIZE << std::endl;
//   int sectionSize = section.first;
//   int sectionNum = section.second;
//   int *trustedM = (int *)malloc(sizeof(int) * M);
//   int j = 0, k;
//   for (int i = 0; i < sectionNum; ++i)
//   {
//     opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId2, 0, 0);
//     k = moveDummy(trustedM, sectionSize);
//     quickSort(trustedM, 0, k - 1);
//     opOneLinearScanBlock(j, trustedM, k, outStructureId1, 1, 0);
//     j += k;
//     if (j > inSize)
//     {
//       std::cout << "Final error2\n";
//     }
//   }
//   std::cout << "Till Final IOcost: " << IOcost / N * BLOCK_DATA_SIZE << std::endl;
//   return outStructureId1;
// }

// std::pair<int, int> ObliviousLooseSort(int inStructureId, int inSize, int outStructureId1, int outStructureId2)
// {
//   std::cout << "In ObliviousLooseSort\n";
//   int *trustedM;
//   if (inSize <= M)
//   {
//     trustedM = (int *)malloc(sizeof(int) * M);
//     opOneLinearScanBlock(0, trustedM, inSize, inStructureId, 0, 0);
//     quickSort(trustedM, 0, inSize - 1);
//     opOneLinearScanBlock(0, trustedM, inSize, outStructureId1, 1, 0);
//     free(trustedM);
//     return {outStructureId1, inSize};
//   }
//   std::vector<int> trustedM2;
//   int realNum = Sample(inStructureId, inSize, trustedM2, is_tight);
//   std::pair<int, int> section = OneLevelPartition(inStructureId, inSize, trustedM2, realNum, P, outStructureId1, 0);
//   int sectionSize = section.first;
//   int sectionNum = section.second;
//   int totalLevelSize = sectionNum * sectionSize;
//   int k;
//   freeAllocate(outStructureId2, outStructureId2, totalLevelSize);
//   trustedM = (int *)malloc(sizeof(int) * M);
//   for (int i = 0; i < sectionNum; ++i)
//   {
//     opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId1, 0, 0);
//     k = moveDummy(trustedM, sectionSize);
//     quickSort(trustedM, 0, k - 1);
//     opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId2, 1, 0);
//   }
//   return {outStructureId2, totalLevelSize};
// }

// // TODO: Finish
// std::pair<int, int> ObliviousLooseSort2(int inStructureId, int inSize, int sampleId, int sortedSampleId, int outStructureId1, int outStructureId2)
// {
//   std::cout << "In ObliviousLooseSort2 && In SasmpleRec\n";
//   std::vector<std::vector<int>> pivots;
//   SampleRec(inStructureId, sampleId, sortedSampleId, 0, pivots);
//   std::cout << "Till Sample IOcost: " << IOcost / N * BLOCK_DATA_SIZE << std::endl;
//   std::cout << "In TwoLevelPartition\n";
//   std::pair<int, int> section = TwoLevelPartition(inStructureId, pivots, P, outStructureId1, outStructureId2);
//   int sectionSize = section.first;
//   int sectionNum = section.second;
//   int totalLevelSize = sectionNum * sectionSize;
//   int k;
//   freeAllocate(outStructureId1, outStructureId1, totalLevelSize);
//   int *trustedM = (int *)malloc(sizeof(int) * M);
//   for (int i = 0; i < sectionNum; ++i)
//   {
//     opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId2, 0, 0);
//     k = moveDummy(trustedM, sectionSize);
//     quickSort(trustedM, 0, k - 1);
//     opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId2, 1, 0);
//   }
//   std::cout << "Till Final IOcost: " << IOcost / N * BLOCK_DATA_SIZE << std::endl;
//   return {outStructureId1, totalLevelSize};
// }

// // TODO: Finish, what's return value
// void ObliviousLooseSortRec(int sampleId, int sampleSize, int sortedSampleId, std::vector<std::vector<int>> &pivots)
// {
//   std::cout << "In ObliviousLooseSortRec\n";
//   std::vector<int> trustedM2;
//   int realNum = Sample(sampleId, sampleSize, trustedM2, 0, 1);
//   std::cout << "In OneLevelPartition\n";
//   std::pair<int, int> section = OneLevelPartition(sampleId, sampleSize, trustedM2, realNum, _P, sortedSampleId, 1);
//   int sectionSize = section.first;
//   int sectionNum = section.second;
//   int j = 0, k = 0, total = 0;
//   int outj = 0, inj = 0;
//   int *trustedM = (int *)malloc(sizeof(int) * M);
//   std::vector<int> quantileIdx;
//   for (int i = 1; i < P; ++i)
//   {
//     quantileIdx.push_back(i * sampleSize / P);
//   }
//   int size = ceil(1.0 * sampleSize / P);
//   std::vector<std::vector<int>> quantileIdx2;
//   std::vector<int> index;
//   for (int i = 0; i < P; ++i)
//   {
//     for (int j = 1; j < P; ++j)
//     {
//       index.push_back(i * size + j * size / P);
//     }
//     quantileIdx2.push_back(index);
//     index.clear();
//   }
//   std::vector<int> pivots1, pivots2_part;
//   for (int i = 0; i < sectionNum; ++i)
//   {
//     opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, sortedSampleId, 0, 0);
//     k = moveDummy(trustedM, sectionSize);
//     quickSort(trustedM, 0, k - 1);
//     total += k;
//     // Cal Level1 pivots
//     while ((j < P - 1) && (quantileIdx[j] < total))
//     {
//       pivots1.push_back(trustedM[quantileIdx[j] - (total - k)]);
//       j += 1;
//     }
//     // Cal Level2 pivots
//     while (outj < P)
//     {
//       while ((inj < P - 1) && (quantileIdx2[outj][inj] < total))
//       {
//         pivots2_part.push_back(trustedM[quantileIdx2[outj][inj] - (total - k)]);
//         inj += 1;
//         if (inj == P - 1)
//         {
//           inj = 0;
//           outj += 1;
//           pivots.push_back(pivots2_part);
//           pivots2_part.clear();
//           break;
//         }
//       }
//       if (outj == P || quantileIdx2[outj][inj] >= total)
//       {
//         break;
//       }
//     }
//     if ((j >= P - 1) && (outj >= P))
//     {
//       break;
//     }
//   }
//   pivots1.insert(pivots1.begin(), INT_MIN);
//   pivots1.push_back(INT_MAX);
//   for (int i = 0; i < P; ++i)
//   {
//     pivots[i].insert(pivots[i].begin(), INT_MIN);
//     pivots[i].push_back(INT_MAX);
//   }
//   pivots.insert(pivots.begin(), pivots1);
// }