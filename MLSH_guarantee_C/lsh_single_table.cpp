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

const int vectorLength = 128;
const int nPlanes = 30;

int main(int argc, char** argv)
{
    int nBaseVecs = 0;    // number of base vectors
    int nQueryVecs = 0;   // number of query vectors


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
    free(hyperplanes);
    free(baseVecs);
    free(queryVecs);
    free(groundTruth);
}