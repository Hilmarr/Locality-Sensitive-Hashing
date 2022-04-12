#include <cmath>
#include <cstring>
#include <iostream>
#include <iterator>
#include <numeric>

#define TIME_LSH
// #define NVTX_PROFILE

#ifdef TIME_LSH
#include <sys/time.h>
#endif

#ifdef NVTX_PROFILE
#include <nvToolsExt.h>
#endif

#include "read_file_funcs.h"

// Globally defined so that the compiler might make assumptions about it
// later on if expedient
const int vectorLength = 128;
const int _THRESHOLD = 18;
const int THRESHOLD = _THRESHOLD * _THRESHOLD;  // Threshold to check on other side of hyperplane(s)

/**
 * Calculates the hash values of the input points given its input hyperplanes,
 * then stores those hash values into an array called indexGroupMap
 * which will later map the point index into a group index used for indexing
 * a second array.
 *
 * @param nPoints : Number of points to hash
 * @param points : The points to hash
 * @param nPlanes : Number of hyperplanes used
 * @param hyperplanes : Array of hyperplanes of length nPlanes
 *
 * @param indexGroupMap : Output - Calculated hash values, may be used as indices
 */
void calculate_hash_values(
    // input
    int nPoints, float* __restrict__ points,
    int nPlanes, float* __restrict__ hyperplanes,
    // output
    int* __restrict__ indexGroupMap)
{
    // - Calculate hash values, keep track of sizes of each group -
#pragma acc data \
  pcopyin(points[nPoints*vectorLength]) \
  pcopyin(hyperplanes[nPlanes*vectorLength]) \
  pcopyout(indexGroupMap[nPoints])
 {
    #pragma acc parallel loop gang worker num_workers(128) vector_length(8)
    for (int i = 0; i < nPoints; i++) {
        float* point = &points[i * vectorLength];
        int hashcode = 0;  //  hashcode will be the group index

        // calculate hash value of the i'th point, store resut in indexGroupMap
        #pragma acc loop reduction(|:hashcode) vector
        for (int j = 0; j < nPlanes; j++) {
            float* hplane = &hyperplanes[j * vectorLength];  // first hyperplane
            // calculate point * hplane
            float vecMul = 0;
            #pragma acc loop reduction(+:vecMul) seq
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
}

/**
 * Calculates the hash values of the input points given its input hyperplanes,
 * then stores those hash values into an array called indexGroupMap
 * which will later map the point index into a group index used for indexing
 * a second array. Also stores the squared distances from each point to each
 * hyperplane in the sqrdDists array
 *
 * @param nPoints : Number of points to hash
 * @param points : The points to hash
 * @param nPlanes : Number of hyperplanes used
 * @param hyperplanes : Array of hyperplanes of length nPlanes
 *
 * @param indexGroupMap : Output - Calculated hash values, may be used as indices
 * @param sqrdDists : Output - Distance from each point to each hyperplane
 */
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
            sqrdDists[i * nPlanes + j] = vecMul * vecMul;
            // printf("vecMul*vecMul = %f\n", vecMul*vecMul);
            // printf("sqrtdDists[%d][%d] = %f\n", i, j, sqrdDists[i * nPlanes + j]);
        }
        indexGroupMap[i] = hashcode;  // save the hashcode
    }
}

/**
 * Takes as input an array of point-index-to-group mappings. Turns this into a
 * data structure which takes the hash value of a point and returns a group index
 * into a new array called groupArray, which itself contains indices into the
 * original point array.
 * The arrays, groupSizeMap and groupIndexMap, are used to easily index or
 * iterate through specific groups in the groupArray. This is handy when we want
 * to compare a point to other points in its group.
 *
 * @param nPoints :            The size of indexGroupMap, groupArray,
 *                             and in general the number of points

 * @param nGroups :            The number of groups.
 *
 * @param indexGroupMap :      Mapping from indices to groups, calculated by function
 *                             'calculate_hash_values'
 *
 * @param groupIndexMapTails : Temporary array used when creating groupArray.
 *                             Preallocated only for possible increased efficiency.
 *
 * @param groupArray    : Output - An array with point indices sorted by their group.
 *                                 Maps into the original point array.
 *
 * @param groupIndexMap : Output - Array of group indices, indexes groupArray
 */
void organize_points_into_groups(
    // input
    int nPoints, int nGroups,
    int* __restrict__ indexGroupMap,
    // output
    int* __restrict__ groupArray,
    int* __restrict__ groupIndexMap)
{
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

// Have to reserve nGroups elements for both combMasks and combDists
// for the function to be completely secure.
// nGroups is the maximum amount of boxes that a point can be mapped to.

/**
 * Takes squared distances as arguments and finds all non-ordered combinations
 * of those distances that are still lower than the threshold. Returns the masks
 * representing the combinations in the combMasks array, and the accumulated
 * distances in the combDists array.
 * 
 * Have to reserve nGroups elements for both combMasks and combDists
 * for the function to be completely secure.
 * nGroups is the maximum amount of boxes that a point can be mapped to.
 *
 * @param sqrdDists : The squared dist for each hyperplane
 * 
 * @param nPlanes : The number of hyperplanes (length of sqrdDists)
 * 
 * @param threshold : The maximum distance allowed squared
 *
 * @param combMasks : Output - Masks where the sum of distances for the sides
 *                             marked as 1 are lower than the threshold.
 *
 * @param combDists : Output - The accumulated sum of distances for each mask
 *
 * @return The length of combMasks
 */
int findNeighborMasks(float* sqrdDists, int nPlanes, float threshold,
                      int* combMasks, float* combDists) {
    // Add single values below the threshold to lists
    int setLen = 0;
    for (int i = 0; i < nPlanes; i++) {
        if (sqrdDists[i] < threshold) {
            combMasks[setLen] = 1 << i;
            combDists[setLen] = sqrdDists[i];
            setLen++;
        }
    }

    // For every value starting with the initial setLen, try to combine the value
    // with all other values that comes after it in the list
    // can be combined if combDists[k] + combDists[i] < threshold
    int underThreshLen = setLen;
    for (int k = underThreshLen - 2; k >= 0; k--) {
        int tmp = setLen;
        for (int i = k + 1; i < setLen; i++) {
            if (combDists[k] < threshold) {
                int dist = combDists[k] + combDists[i];
                if (dist < threshold) {
                    combMasks[tmp] = combMasks[k] | combMasks[i];
                    combDists[tmp] = dist;
                    tmp++;
                }
            }
        }
        setLen = tmp;
    }

    return setLen;
}

/**
 * For every point it goes through all the groups in groupArray it can get mapped into.
 * groupArray stores these points as indices into the original point array,
 * it does this for all tables.
 * Then the function puts indices it can get mapped into that point's place in
 * the potentialMatches array. After that the function creates two other arrays
 * called potentialMatchesIndices and potentialMatchesLengths; these two arrays
 * are used to index and iterate the potentialMatches array.
 *
 * @param nPoints: Number of points in the LSH tables and size of groupArray
 * @param indexGroupMap : Calculated hash values
 * @param indexGroupMapTableLen : Length of indexGroupMap per table
 * @param groupArray : Arrays with point indices sorted by their group.
 * @param groupIndexMap : Arrays with group indices, indexes into groupArray.
 *
 * @param pPotentialMatches : Output.
 *              Pointer to an array with indices into point1 containing the
 *              indices of all possible matches for each point in points2
 *
 * @param pPotentialMatchesIndices : Output.
 *              Pointer to array with indices into potentialMatches
 *
 * @return The total amount of matches
 */
int find_potential_matches(
    int nPlanes,
    int nGroups,
    int nPoints,
    int* __restrict__ indexGroupMap,
    float* __restrict__ sqrdDists,
    int* __restrict__ groupArray,
    int* __restrict__ groupIndexMap,
    int** __restrict__ pPotentialMatches,
    int** __restrict__ pPotentialMatchesIndices)
{
    // int potentialMatchesLen = 1e4 * nPoints;
    int potentialMatchesLen = 1e9;
    int* potentialMatches = (int*)malloc(potentialMatchesLen * sizeof(int));
    int* potentialMatchesIndices = (int*)malloc((nPoints + 1) * sizeof(int));

    // Masks and distances for all combinations of hyperplanes where the euclidean
    // distance from the point to the intersection of the hyperplane is below
    // a set threshold.
    // nGroups is the theoretical upper limit of possible combinations, but
    // in reality the number of potential matches will likely never be that high
    int* combMasks = (int*)malloc((1+nGroups) * sizeof(int));
    combMasks[0] = 0;
    float* combDists = (float*)malloc(nGroups * sizeof(float));
    int cnt = 0;
    for (int i = 0; i < nPoints; i++) {
        // Set the index to the groups that the query point is mapped to
        potentialMatchesIndices[i] = cnt;
        // Get the hashcode for this query point
        int hashcode = indexGroupMap[i];
        // Find all combinations of hyperplanes that is closer than sqrt(THRESHOLD)
        int setLen = findNeighborMasks(&sqrdDists[i * nPlanes], nPlanes, THRESHOLD,
                                       combMasks+1, combDists);
        // Find potential matches by checking hashcode and all hashcodes within
        // a certain distance
        for (int j = 0; j < setLen+1; j++) {
            // get next hashcode
            int hc = hashcode ^ combMasks[j];
            // Find start and end index of the group
            int kStart = groupIndexMap[hc];
            int kEnd = groupIndexMap[hc+1];
            int groupSize = kEnd - kStart;
            // Check if we have reserved enough memory
            if (cnt + groupSize >= potentialMatchesLen) {
                fprintf(stderr, "Not enough with %d elements for potentialMatches. ", potentialMatchesLen);
                exit(1);
            }
            // Add point indices in group hc to potentialMatches
            for (int k = kStart; k < kEnd; k++) {
                potentialMatches[cnt] = groupArray[k];
                cnt++;
            }
        }
    }
    potentialMatchesIndices[nPoints] = cnt;

    free(combMasks);
    free(combDists);

    *pPotentialMatches = potentialMatches;
    *pPotentialMatchesIndices = potentialMatchesIndices;

    return cnt;
}

/**
 * Uses result from find_potential_matches to match points1 and points2
 *
 * @param nQueryVecs : Number of query vectors
 *
 * @param nBaseVecs : Number of base vectors
 *
 * @param queryVecs : Query vectors
 *
 * @param baseVecs : Base vectors (Array of points used to create the LSH tables)
 * 
 * @param nPotentialMatches : Number of potential matches
 *
 * @param potentialMatches : An array with indices into point1 containing the
 *                           indices of all possible matches for each point in points2
 *                           (i.e. potential matches between points2 and points1)
 *
 * @param potentialMatchesIndices : Indices used to index into potentialMatches
 *
 *
 * @param lshMatches : Output - Index into points1 representing the best match
 *                              between points2  and points1
 *
 * @param bestMatchDists : Output - Distances between points for every match
 *                                  in lshMatches
 *
 * @param lshMatches2 : Output - Index into points1 representing the 2nd best match
 *                              between points2  and points1
 *
 * @param bestMatchDists2 : Output - Distances between points for every match
 *                                  in lshMatches2
 */
void match_points(int nQueryVecs,
                  int nBaseVecs,
                  float* __restrict__ queryVecs,
                  float* __restrict__ baseVecs,
                  int nPotentialMatches,
                  int* __restrict__ potentialMatches,
                  int* __restrict__ potentialMatchesIndices,
                  int* __restrict__ potentialMatchesLengths,
                  // outputs
                  int* __restrict__ lshMatches, float* __restrict__ bestMatchDists,
                  int* __restrict__ lshMatches2, float* __restrict__ bestMatchDists2)
{
#pragma acc data \
    pcopyin(queryVecs[nQueryVecs*vectorLength]) \
    pcopyin(baseVecs[nBaseVecs*vectorLength]) \
    pcopyin(potentialMatches[nPotentialMatches])\
    pcopyin(potentialMatchesIndices[nQueryVecs]) \
    pcopyin(potentialMatchesLengths[nQueryVecs]) \
    pcopy(lshMatches[nQueryVecs]) \
    pcopy(bestMatchDists[nQueryVecs]) \
    pcopy(lshMatches2[nQueryVecs]) \
    pcopy(bestMatchDists2[nQueryVecs]) 
 {
    #pragma acc parallel loop  num_workers(1) vector_length(32)
    for (int i = 0; i < nQueryVecs; i++) {
        // Find the group of elements from groupArray to match with
        bool changed = false;
        float bestMatchDist = bestMatchDists[i];
        float bestMatchDist2 = bestMatchDist;
        int match = -1;
        int match2 = -1;
        int jStart = potentialMatchesIndices[i];
        int jEnd = jStart + potentialMatchesLengths[i];

        // Match the points
        for (int j = jStart; j < jEnd; j++) {
            // printf("u=%d\n", u);
            int idx = potentialMatches[j];
            
            // diff = sum((queryVecs[i][:] - baseVecs[idx][:]) .^ 2);
            float diff = 0;
            #pragma acc loop reduction(+:diff) vector
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
                changed = true;
            }
        }

        if (changed)
        {
            lshMatches[i] = match;
            bestMatchDists[i] = bestMatchDist;
            lshMatches2[i] = match2;
            bestMatchDists2[i] = bestMatchDist2;
        }
    }
 }
}


int main(int argc, char** argv) {
    int nBaseVecs = 0;   // number of base vectors
    int nQueryVecs = 0;  // number of query vectors

    const int nPlanes = 16;

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

#ifdef TIME_LSH
    struct timeval time;
    long startTime;
    long endTime;
    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

#ifdef NVTX_PROFILE
    nvtxRangePush("Calculating hash values for base vectors");
#endif

    calculate_hash_values(nBaseVecs, baseVecs, nPlanes, hyperplanes,
                          indexGroupMap);

#ifdef NVTX_PROFILE
    nvtxRangePop();
#endif

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    printf("calculate_hash_values for base vectors:\n");
    printf("   - time: %.3f seconds\n", ((double)endTime - startTime) / 1000);

    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

#ifdef NVTX_PROFILE
    nvtxRangePush("Constructing lsh tables");
#endif

    organize_points_into_groups(nBaseVecs, nGroups, indexGroupMap,
                                groupArray, groupIndexMap);

#ifdef NVTX_PROFILE
    nvtxRangePop();
#endif

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    printf("organize_points_into_groups:\n");
    printf("   - time: %.3f seconds\n", ((double)endTime - startTime) / 1000);

    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

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

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

#ifdef NVTX_PROFILE
    nvtxRangePush("Calculating hash values for query vectors");
#endif

    calculate_hash_values_and_dists(nQueryVecs, queryVecs, nPlanes, hyperplanes,
                                    indexGroupMap, sqrdDists);

#ifdef NVTX_PROFILE
    nvtxRangePop();
#endif

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    printf("calculate_hash_values_and_dists for query vectors:\n");
    printf("   - time: %.3f seconds\n", ((double)endTime - startTime) / 1000);

    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

    int* potentialMatches;
    int* potentialMatchesIndices;

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

#ifdef NVTX_PROFILE
    nvtxRangePush("Finding potential matches");
#endif

    int nPotentialMatches;

    nPotentialMatches = find_potential_matches(
        nPlanes, nGroups, nQueryVecs,
        indexGroupMap, sqrdDists,
        groupArray, groupIndexMap,
        &potentialMatches, &potentialMatchesIndices);

#ifdef NVTX_PROFILE
    nvtxRangePop();
#endif

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    printf("find_potential_matches for query vectors:\n");
    printf("   - time: %.3f seconds\n", ((double)endTime - startTime) / 1000);

    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

    free(indexGroupMap);

    int* potentialMatchesLengths = (int*) malloc(nQueryVecs * sizeof(int));
    for (int i = 0; i < nQueryVecs; i++) {
        potentialMatchesLengths[i] = potentialMatchesIndices[i+1] - potentialMatchesIndices[i];
    }


#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

#ifdef NVTX_PROFILE
    nvtxRangePush("Matching potential matches");
#endif

    // - Compare points -
    match_points(nQueryVecs, nBaseVecs, queryVecs, baseVecs, nPotentialMatches,
                 potentialMatches, potentialMatchesIndices, potentialMatchesLengths,
                 lshMatches, bestMatchDists, lshMatches2, bestMatchDists2);

#ifdef NVTX_PROFILE
    nvtxRangePop();
#endif

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    printf("match points:\n");
    printf("   - time: %.3f seconds\n", ((double)endTime - startTime) / 1000);

    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

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

    printf("\n");
    printf("Potential matches found: %d\n", nPotentialMatches);
    double matchesPerQueryVector = ((double)nPotentialMatches) / nQueryVecs;
    printf("Comparisons per query vector: %f\n", matchesPerQueryVector);
    printf("Average portion of search space searched: %f\n", matchesPerQueryVector / nBaseVecs);
    printf("\n");

    // --- Free memory ---

    free(sqrdDists);
    free(groupIndexMap);
    free(groupArray);
    free(potentialMatches);
    free(potentialMatchesIndices);
    free(potentialMatchesLengths);

    free(lshMatches);
    free(bestMatchDists);
    free(lshMatches2);
    free(bestMatchDists2);

    free(hyperplanes);
    free(baseVecs);
    free(queryVecs);
    free(groundTruth);
}
