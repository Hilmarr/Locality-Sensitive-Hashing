#include <cmath>
#include <cstring>
#include <iostream>
#include <iterator>
#include <numeric>

#define TIME_LSH

#ifdef TIME_LSH
#include <sys/time.h>
#endif

#include "read_vector_file_funcs.h"

// Globally defined so that the compiler might make assumptions about it
// later on if expedient
const int vectorLength = 128;

void calculate_hash_values(
    // input
    int nPoints, float* __restrict__ points,
    int nPlanes, float* __restrict__ hyperplanes,
    // output
    int* __restrict__ indexGroupMap)
{
    // - Calculate hash values, keep track of sizes of each group -
    for (int i = 0; i < nPoints; i++) {

        float* point = &points[i * vectorLength];
        int hashcode = 0;  //  hashcode will be the group index

        // calculate hash value of the i'th point, store resut in indexGroupMap
        for (int j = 0; j < nPlanes; j++) {
            float* hplane = &hyperplanes[j * vectorLength];  // first hyperplane
            // calculate point * hplane
            float vecMul = 0;
            for (int k = 0; k < vectorLength; k++) {
                vecMul += point[k] * hplane[k];
            }
            // set i'th bit to one if point is on "positive" side of hyperplane
            if (vecMul > 0) {
                hashcode = hashcode | (1 << j);
            }
        }
        indexGroupMap[i] = hashcode;  // save the hashcode
    }
}

void organize_points_into_groups(
    // input
    int nPoints, int nGroups,
    int* __restrict__ indexGroupMap,
    // output
    int* __restrict__ groupArray,
    int* __restrict__ groupIndexMap)
{
    int* groupSizeMap = (int*) calloc(nGroups, sizeof(int));
    int* groupIndexMapTails = (int*) calloc(nGroups, sizeof(int));

    // find the size of each group
    for (int i = 0; i < nPoints; i++) {
        groupSizeMap[indexGroupMap[i]]++;
    }

    // Find group indices into the group array 'groupArray'
    // // find group indices using exclusive scan of the group sizes
    //  std::exclusive_scan(groupSizeMap, groupSizeMap + nGroups - 1);
    // doing it manually for now
    int cnt = 0;
    for (int i = 0; i < nGroups; i++) {
        groupIndexMap[i] = cnt;
        groupIndexMapTails[i] = cnt;
        cnt += groupSizeMap[i];
        // printf("cnt = %d\n", cnt);
    }
    groupIndexMap[nGroups] = cnt;

    // Fill groupArray
    for (int i = 0; i < nPoints; i++) {
        // what group did the point get hashed into
        int hashcode = indexGroupMap[i];
        // start index of that group in groupArray + number of elements currently inserted
        int idx = groupIndexMapTails[hashcode]++;
        // add the point to the group
        groupArray[idx] = i;
    }
    
    free(groupSizeMap);
    free(groupIndexMapTails);
}

int main(int argc, char** argv)
{
    int nBaseVecs = 0;    // number of base vectors
    int nQueryVecs = 0;   // number of query vectors

    const int nPlanes = 8;

    // -----    Read data    -----

    // Check number of arguments
    if (argc < 5) {
        fprintf(stderr, "Program needs to be given 3 files as arguments:\n");
        fprintf(stderr, " 1. .fvec vector file with feature descriptors (base set)\n");
        fprintf(stderr, " 2. .fvec vector file with feature descriptors (query set)\n");
        fprintf(stderr, " 3. .ivec vector file (Ground truth)\n");
        fprintf(stderr, " 4. binary file with numbers for the hyperplanes\n");
        return -1;
    }

    // - Read hyperplanes (assume there is enough for now) -
    // hyperplanes should be orthogonal otherwise the algorithm won't work
    // properly
    float* hyperplanes = (float*)malloc(nPlanes * vectorLength * sizeof(float));
    FILE* fp = fopen(argv[4], "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error opening %s for reading\n", argv[4]);
        return -1;
    }
    int ret = fread(hyperplanes, 4, nPlanes * vectorLength, fp);
    if (ret != nPlanes * vectorLength) {
        printf("Only %d of %d elements read\n", ret, nPlanes * vectorLength);
        fclose(fp);
        free(hyperplanes);
        return 1;
    }
    fclose(fp);

    // - Read base and query vectors -
    int* tmp1;
    int* tmp2;
    nBaseVecs = read_vector_file2(argv[1], &tmp1) / 128;
    nQueryVecs = read_vector_file2(argv[2], &tmp2) / 128;
    float* baseVecs = (float*)tmp1;
    float* queryVecs = (float*)tmp2;

    // TODO: make LSH tables

    //  -- Arrays to organize groups --
    const int nGroups = 1 << nPlanes;
    // the group that each point falls into
    int* indexGroupMap = (int*)malloc(nBaseVecs * sizeof(int));
    // // distance from each hyperplane for each point
    // float* distanceFromHP = (float*) malloc(nBaseVecs * nPlanes * sizeof(float));

    int* groupIndexMap = (int*) malloc((nGroups + 1) * sizeof(int));
    int* groupArray = (int*) malloc(nBaseVecs * sizeof(int));

    calculate_hash_values(nBaseVecs, baseVecs, nPlanes, hyperplanes,
                          indexGroupMap);

    organize_points_into_groups(nBaseVecs, nGroups, indexGroupMap,
                                groupArray, groupIndexMap);

    // TODO: match

    // - Read ground truth -
    int* tmpGT;
    int groundTruthLen = read_vector_file2(argv[3], &tmpGT) / 100;
    int* groundTruth = (int*)malloc(groundTruthLen * sizeof(int));
    for (int i = 0; i < groundTruthLen; i++) {
        groundTruth[i] = tmpGT[100 * i];
    }
    free(tmpGT);

    // TODO: check if matches are correct


    // free memory
    
    free(indexGroupMap);
    // free(distanceFromHP);
    free(groupIndexMap);
    free(groupArray);

    free(hyperplanes);
    free(baseVecs);
    free(queryVecs);
    free(groundTruth);
}