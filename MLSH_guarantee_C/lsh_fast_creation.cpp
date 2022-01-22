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
const int _THRESHOLD = 20;
const int THRESHOLD = _THRESHOLD * _THRESHOLD; // Threshold to check on other side of hyperplane(s)


inline int set_int_arr_size(int** arr, int curSize, int newSize) {
    int* oldArr = *arr;
    *arr = (int*)malloc(newSize * sizeof(int));
    memcpy(*arr, oldArr, curSize * sizeof(int));
    free(oldArr);
    return newSize;
}

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

void calculate_hash_values_and_dists(
    // input
    int nPoints, float* __restrict__ points,
    int nPlanes, float* __restrict__ hyperplanes,
    // output
    int* __restrict__ indexGroupMap,
    float* __restrict__ sqrdDists)
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
            sqrdDists[i * nPlanes + j] = vecMul;
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
    int* __restrict__ groupIndexMap) {
    int* groupSizeMap = (int*)calloc(nGroups, sizeof(int));
    int* groupIndexMapTails = (int*)calloc(nGroups, sizeof(int));

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

void match_points(int nQueryVecs,
                  float* __restrict__ queryVecs,
                  float* __restrict__ baseVecs, int* __restrict__ indexGroupMap,
                  int* __restrict__ groupIndexMap, int* __restrict__ groupArray,
                  // outputs
                  int* __restrict__ lshMatches, float* __restrict__ bestMatchDists,
                  int* __restrict__ lshMatches2, float* __restrict__ bestMatchDists2) {
    for (int i = 0; i < nQueryVecs; i++) {
        // Find the group of elements from groupArray to match with
        float bestMatchDist = 1e10;
        float bestMatchDist2 = 1e10;
        int match = -1;
        int match2 = -1;
        int hashcode = indexGroupMap[i];
        int jStart = groupIndexMap[hashcode];
        int jEnd = groupIndexMap[hashcode + 1];

        // Match the points
        for (int j = jStart; j < jEnd; j++) {
            int idx = groupArray[j];
            // diff = sum((queryVecs[i][:] - baseVecs[idx][:]) .^ 2);
            float diff = 0;
            for (int k = 0; k < vectorLength; k++) {
                float tmp = queryVecs[i * vectorLength + k] - baseVecs[idx * vectorLength + k];
                diff += tmp * tmp;
            }
            // check if distance is the lowest distance so far
            if (diff < bestMatchDist) {
                // save old value as second best
                bestMatchDist2 = bestMatchDist;
                match2 = match;
                // update best value
                bestMatchDist = diff;
                match = idx;
            }
        }

        lshMatches[i] = match;
        bestMatchDists[i] = bestMatchDist;
        lshMatches2[i] = match2;
        bestMatchDists2[i] = bestMatchDist2;
    }
}

// Have to reserve nGroups elements for both combMasks and combDists
// for the function to be completely secure.
// nGroups is the maximum amount of boxes that a point can be mapped to.
int findNeighborMasks(float* sqrdDists, int nPlanes, float threshold,
                      int* combMasks, float* combDists)
{
    // Add single values below the threshold to lists
    int setLen = 0;
    for (int i = 0; i < nPlanes; i++) {
        if (sqrdDists[i] < threshold) {
            combMasks[setLen] = 1 << i;
            combDists[setLen] = sqrdDists[i];
        }        
    }
    
    // For every value starting with the initial setLen, try to combine the value
    // with all other values that comes after it in the list
    // can be combined if combDists[k] + combDists[i] < threshold
    int underThreshLen = setLen;
    for (int k = underThreshLen-2; k >= 0; k--) {
        int tmp = setLen;
        for (int i = k+1; i < setLen; i++) {
            if (combDists[k] < threshold) {
                int dist = combDists[k] + combDists[i];
                if (dist < threshold) {
                    tmp++;
                    combMasks[tmp] = combMasks[k] | combMasks[i];
                    combDists[tmp] = dist;
                }
            }
        }
        setLen = tmp;
    }

    return setLen;
}


int main(int argc, char** argv) {
    int nBaseVecs = 0;   // number of base vectors
    int nQueryVecs = 0;  // number of query vectors

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

    // ---  Make LSH tables ---

    //  - Arrays to organize groups -
    const int nGroups = 1 << nPlanes;
    // the group that each base point falls into
    int* indexGroupMap = (int*)malloc(nBaseVecs * sizeof(int));

    int* groupIndexMap = (int*)malloc((nGroups + 1) * sizeof(int));
    int* groupArray = (int*)malloc(nBaseVecs * sizeof(int));

    calculate_hash_values(nBaseVecs, baseVecs, nPlanes, hyperplanes,
                          indexGroupMap);

    organize_points_into_groups(nBaseVecs, nGroups, indexGroupMap,
                                groupArray, groupIndexMap);

    free(indexGroupMap);

    // --- Match points ---

    // - Array allocation and initialization -
    // holds the actual matches
    int* lshMatches = (int*)malloc(nQueryVecs * sizeof(int));
    memset(lshMatches, -1, nQueryVecs * sizeof(int));
    // holds the best match distance for each query vector
    float* bestMatchDists = (float*)malloc(nQueryVecs * sizeof(float));
    for (int i = 0; i < nQueryVecs; i++) {
        bestMatchDists[i] = 1e10;
    }

    // Holds the second best matches
    int* lshMatches2 = (int*)malloc(nQueryVecs * sizeof(int));
    memset(lshMatches2, -1, nQueryVecs * sizeof(int));
    // holds the second best match distance for each query vector
    float* bestMatchDists2 = (float*)malloc(nQueryVecs * sizeof(float));
    for (int i = 0; i < nQueryVecs; i++) {
        bestMatchDists[i] = 1e10;
    }

    // squared distance from each hyperplane for each point
    float* sqrdDists = (float*)malloc(nBaseVecs * nPlanes * sizeof(float));

    // the group that each query point falls into
    indexGroupMap = (int*)malloc(nQueryVecs * sizeof(int));

    // - Calculate hash points for query vectors -

    calculate_hash_values_and_dists(nQueryVecs, queryVecs, nPlanes, hyperplanes,
                          indexGroupMap, sqrdDists);

    int idxGroupMapExtLen = 5;
    int* idxGroupMapExt = (int*) malloc(idxGroupMapExtLen * sizeof(int));
    int* idxGroupMapExtIndices = (int*)malloc((nQueryVecs+1) * sizeof(int));

    // Masks and distances for all combinations of hyperplanes where the euclidean
    // distance from the point to the intersection of the hyperplane is below
    // a set threshold.
    // nGroups is the theoretical upper limit of possible combinations, but
    // in reality the number of potential matches will likely never be that high
    int* combMasks = (int*) malloc(nGroups * sizeof(int));
    float* combDists = (float*) malloc(nGroups * sizeof(float));

    int cnt = 0;
    for (int i = 0; i < nQueryVecs; i++) {
        // Set the index to the groups that the query point is mapped to
        idxGroupMapExtIndices[i] = cnt;
        // Get the hashcode for this query point
        int hashcode = indexGroupMap[i];
        // Add hashcode to the list
        idxGroupMapExt[cnt] = hashcode;
        cnt++;
        // Find all combinations of hyperplanes that is closer than sqrt(THRESHOLD)
        int setLen = findNeighborMasks(sqrdDists, nPlanes, THRESHOLD,
                                       combMasks, combDists);
        // Check if we have reserved enough memory
        if (cnt + setLen >= idxGroupMapExtLen) {
            printf("Not enough with %d elements for idxGroupMapExt. ", idxGroupMapExtLen);
            printf("Expanding to %d elements\n", 2*(cnt+setLen) + 1);
            idxGroupMapExtLen = set_int_arr_size(
                &idxGroupMapExt, idxGroupMapExtLen, (cnt+setLen) * 2 + 1);
        }
        // Add neighboring boxes to the extended groupMap
        for (int j = 0; j < setLen; j++) {
            idxGroupMapExt[cnt] = hashcode ^ combMasks[j];
            cnt++;
        }
    }
    idxGroupMapExtIndices[nQueryVecs] = cnt;

    free(combMasks);
    free(combDists);

    // - Compare points -
    match_points(nQueryVecs, queryVecs, baseVecs, indexGroupMap,
                 groupIndexMap, groupArray,
                 lshMatches, bestMatchDists, lshMatches2, bestMatchDists2);

    free(indexGroupMap);

    // --- Check matches ---

    // - Read ground truth -
    int* tmpGT;
    int groundTruthLen = read_vector_file2(argv[3], &tmpGT) / 100;
    int* groundTruth = (int*)malloc(groundTruthLen * sizeof(int));
    for (int i = 0; i < groundTruthLen; i++) {
        groundTruth[i] = tmpGT[100 * i];
    }
    free(tmpGT);

    // - Check how many matches were correct -
    int correct = 0;
    for (int i = 0; i < nQueryVecs; i++) {
        correct += lshMatches[i] == groundTruth[i];
    }
    double correctRatio = ((double)correct) / nQueryVecs;
    printf("Correct ratio: %f\n", correctRatio);

    long unsigned int diff_both = 0;
    long unsigned int diff_correct = 0;
    long unsigned int diff_incorrect = 0;
    for (int i = 0; i < nQueryVecs; i++) {
        double diff = 0;
        float* p_query = &queryVecs[i * 128];
        float* p_base = &baseVecs[groundTruth[i] * 128];
        for (int j = 0; j < 128; j++) {
            diff += (p_query[i] - p_base[i]) * (p_query[i] - p_base[i]);
        }
        diff = sqrt(diff);
        diff_both += diff;
        if (lshMatches[i] == groundTruth[i]) {
            diff_correct += (int)diff;
        } else {
            diff_incorrect += (int)diff;
        }
    }
    int incorrect = nQueryVecs - correct;

    printf("Average distance of best fits: %f\n",
           ((double)diff_both) / nQueryVecs);
    printf("Average distance correctly classified points: %f\n",
           ((double)diff_correct) / correct);
    printf("Average distance incorrectly classified points: %f\n",
           ((double)diff_incorrect) / incorrect);

    // --- Free memory ---

    free(sqrdDists);
    free(groupIndexMap);
    free(groupArray);
    free(idxGroupMapExt);
    free(idxGroupMapExtIndices);

    free(lshMatches);
    free(bestMatchDists);
    free(lshMatches2);
    free(bestMatchDists2);

    free(hyperplanes);
    free(baseVecs);
    free(queryVecs);
    free(groundTruth);
}
