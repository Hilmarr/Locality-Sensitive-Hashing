#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <cstring>

#define TIME_LSH

#ifdef TIME_LSH
#include <sys/time.h>
#endif

/**
 * Returns a random float between -1 and 1
 */
inline float rand_minus1_to_1() {
    return ((float)rand() / (RAND_MAX/2)) - 1;
}

/**
 * Fill points1 with random values.
 * Fills points2 with values in points1 + some uniform noise.
 * 
 * @param nPoints1 : Length of points1 array
 * @param nPoints2 : Length of points2 array
 * @param vectorLength : Number of dimensions for each point
 * @param noisScale : The amount of noise relative to the length of each dimension
 * 
 * @param points1 : Output - Point array 1
 * @param points2 : Output - Point array 2
 */
void fill_point_arrays(int nPoints1, int nPoints2, int vectorLength, float noiseScale,
                       float* __restrict__ points1, float* __restrict__ points2);

/**
 * Fills the hyperplane array with random hyperplane vectors.
 * 
 * @param nPlanes : Number of hyperplanes
 * @param vectorLength : Number of dimensions for each hyperplane
 * @param hyperplanes: Output - The returned hyperplanes
 */
void fill_hyperplanes(int nPlanes, int vectorLength, float* hyperplanes);

/**
 * Calculates the hash values of the input points given its input hyperplanes,
 * then stores those hash values into an array called indexGroupMap
 * which will later map the point index into a group index used for indexing
 * a second array.
 * 
 * @param vectorLength : Number of dimensions for each point
 * @param nPoints : Number of points to hash
 * @param points : The points to hash
 * @param nPlanes : Number of hyperplanes used
 * @param hyperplanes : Array of hyperplanes of length nPlanes
 * 
 * @param indexGroupMap : Output - Calculated hash values, may be used as indices
 */
void calculate_hash_values(
    // input
    int vectorLength,
    int nPoints, float* __restrict__ points,
    int nPlanes, float* __restrict__ hyperplanes,
    // output
    int* __restrict__ indexGroupMap);

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
 * @param groupSizeMap  : Output - Array of group sizes
 * 
 * @param groupIndexMap : Output - Array of group indices, indexes groupArray
 */
void organize_points_into_groups(
    // input
    int nPoints, int nGroups,
    int* __restrict__ indexGroupMap,
    // preallocated temporary storage
    int* __restrict__ groupIndexMapTails,
    // output
    int* __restrict__ groupArray,
    int* __restrict__ groupSizeMap,
    int* __restrict__ groupIndexMap);

/**
 * Calculates indexGroup for all tables by calling 'calculate_hash_values'
 * int a loop.
 * 
 * @param vectorLength : Number of dimensions for each point
 * @param numTables : Number of LSH tables
 * @param nPoints : Number of points to hash
 * @param points : The points to hash
 * @param nPlanes : Number of hyperplanes used
 * @param hyperplanes : Array of hyperplanes of length nPlanes
 * @param indexGroupMapTableLen : The length of a single table in indexGroupMap
 * 
 * @param indexGroupMap : Output - Calculated hash values, may be used as indices
 */
void calculate_indexGroupMap(
    // input
    int vectorLength, int numTables,
    int nPoints, float* __restrict__ points,
    int nPlanes, float* __restrict__ hyperplanes,
    int indexGroupMapTableLen,
    // output
    int* __restrict__ indexGroupMap);

/**
 * Combines functionality of 'calculate_indexGroupMap' and 'organize_points_into_group'
 * in order to create a set of LSH tables which can later be used to match the
 * dataset with other datasets.
 * 
 * @param vectorLength : Number of dimensions for each point
 * @param numTables : Number of LSH tables
 * @param nGroups : Number of groups
 * @param nPoints : Number of points to hash
 * @param points : The points to hash
 * @param nPlanes : Number of hyperplanes used
 * @param hyperplanes : Array of hyperplanes of length nPlanes
 * @param indexGroupMap : Calculated hash values, i.e. indices into groups
 * @param indexGroupMapTableLen : The length of a single table in indexGroupMap
 * 
 * @param groupArray : Output - Arrays with point indices sorted by their group.
 *                              Maps into the original point array.
 * 
 * @param groupSizeMap : Output - Arrays with sizes of each group in groupArray.
 * 
 * @param groupIndexMap : Output - Arrays with group indices, indexes into groupArray.
 * 
 */
void construct_lsh_tables(  // input
    int vectorLength, int numTables, int nGroups,
    int nPoints, float* __restrict__ points,
    int nPlanes, float* __restrict__ hyperplanes,
    int* __restrict__ indexGroupMap, int indexGroupMapTableLen,
    // output
    int* __restrict__ groupArray,
    int* __restrict__ groupSizeMap,
    int* __restrict__ groupIndexMap);

/**
 * For every point it goes through all the groups in groupArray it can get mapped into.
 * groupArray stores these points as indices into the original point array,
 * it does this for all tables.
 * Then the function puts indices it can get mapped into that point's place in
 * the potentialMatches array. After that the function creates two other arrays
 * called potentialMatchesIndices and potentialMatchesLengths; these two arrays
 * are used to index and iterate the potentialMatches array.
 * 
 * @param numTables : Number of LSH tables
 * @param nPoints1 : Number of points in the LSH tables (points1)
 *                   (size of groupArray, groupSizeMap, and groupIndexMap)
 * @param nPoints2 : Number of points to be matched with the points in the LSH tables (points2)
 * @param indexGroupMap : Calculated hash values for points2
 * @param indexGroupMapTableLen : Length of indexGroupMap per table
 * @param groupArray : Arrays with point indices sorted by their group.
 * @param groupSizeMap : Arrays with sizes of each group in groupArray.
 * @param groupIndexMap : Arrays with group indices, indexes into groupArray.
 * @param groupMapTableLen : Length of all the groupMaps (number of groups)
 * @param potentialMatchesMaxLen : The current max size of potentialMatchesLen
 * 
 * @param potentialMatches : Output - An array with indices into point1 containing the
 *                           indices of all possible matches for each point in points2
 * 
 * @param potentialMatchesIndices : Output - Indices into potentialMatches
 * 
 * @param potentialMatchesLengths : Output - Number of matches for each point in point2
 *                                           Used to iterate through potentialMatches
 * 
 * @return The total amount of matches
 */
int find_potential_matches(
    // inputs
    int numTables, int nPoints1, int nPoints2,
    int* __restrict__ indexGroupMap, int indexGroupMapTableLen,
    int* __restrict__ groupArray, int* __restrict__ groupSizeMap,
    int* __restrict__ groupIndexMap,
    int groupMapTableLen, int potentialMatchesMaxLen,
    // outputs
    int** __restrict__ potentialMatches,
    int* __restrict__ potentialMatchesIndices,
    int* __restrict__ potentialMatchesLengths);

/**
 * Uses result from find_potential_matches to match points1 and points2
 * 
 * @param vectorLength : Number of dimensions for each point
 * 
 * @param numTables : Number of LSH tables
 * 
 * @param nPoints2 : Number of points in points 2
 * 
 * @param points2 : Array of points to be matched with points1
 * 
 * @param points1 : Array of points used to create the LSH tables
 * 
 * @param potentialMatches : An array with indices into point1 containing the
 *                           indices of all possible matches for each point in points2
 *                           (i.e. potential matches between points2 and points1)
 *
 * @param potentialMatchesIndices : Indices used to index into potentialMatches
 *
 * @param potentialMatchesLengths : Number of matches for each point in point2
 *                                  Used to iterate through potentialMatches
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
void match_points(
    // inputs
    int vectorLength, int nPoints2,
    float* __restrict__ points2, float* __restrict__ points1,
    int* __restrict__ potentialMatches, int* __restrict__ potentialMatchesIndices,
    int* __restrict__ potentialMatchesLengths,
    // outputs
    int* __restrict__ lshMatches, float* __restrict__ bestMatchDists,
    int* __restrict__ lshMatches2, float* __restrict__ bestMatchDists2);

int double_int_arr_size(int** arr, int curSize) {
    int newSize = 2*curSize;
    int* newArr = (int*) malloc(newSize * sizeof(int));
    for (int i = 0; i < curSize; i++) {
        newArr[i] = (*arr)[i];
    }
    int* oldArr = *arr;
    *arr = newArr;
    free(oldArr);
    return newSize;
}


int main() {
    srand(0);
    const int nPoints1 = 40000;     // number of points in the first dataset
    const int nPoints2 = nPoints1;  // number of points in the second dataset
    const int vectorLength = 128;
    const float noiseScale = 0.3;
    const int numTables = 16;

    // -- Generate some random points and similar random points to match with --

    // points to be matched
    float* points1 = (float*)malloc(numTables * nPoints1 * vectorLength * sizeof(float));
    float* points2 = (float*)malloc(numTables * nPoints2 * vectorLength * sizeof(float));

    const int nPlanes = (int)log2(nPoints1);

    fill_point_arrays(nPoints1, nPoints2, vectorLength, noiseScale, points1, points2);

    // -- Generate hyperplanes --

    // allocate hyperplanes
    const int hyperplanesTableLen = nPlanes * vectorLength;
    const int hyperplanesLen = numTables * hyperplanesTableLen;
    float* hyperplanes = (float*)malloc(hyperplanesLen * sizeof(float));

    // fill all hyperplanes
    for (int table = 0; table < numTables; table++) {
        float* hyperplanes2 = hyperplanes + table * hyperplanesTableLen;
        // create normalized hyperplanes represented by vectors of euclidean size 1
        fill_hyperplanes(nPlanes, vectorLength, hyperplanes2);
    }

    //  -- Arrays to organize groups --

    const int nGroups = 1 << nPlanes;
    // the group that each point falls into
    const int indexGroupMapTableLen = ((nPoints1 > nPoints2) ? nPoints1 : nPoints2) * sizeof(int);
    const int indexGroupMapLen = numTables * indexGroupMapTableLen;
    int* indexGroupMap = (int*)malloc(indexGroupMapLen * sizeof(int));

    // in case we're changing sizes
    const int groupMapTableLen = nGroups;
    const int groupMapLen = numTables * groupMapTableLen;
    // amount of points that fall into each group
    int* groupSizeMap = (int*)malloc(groupMapLen * sizeof(int));
    // the starting index for each group in groupArray
    int* groupIndexMap = (int*)malloc(groupMapLen * sizeof(int));

    // Actual groups
    const int groupArrayTableLen = nPoints1;
    const int groupArrayLen = numTables * groupArrayTableLen;
    int* groupArray = (int*)malloc(groupArrayLen * sizeof(int));


#ifdef TIME_LSH
    struct timeval time;
    long startTime;
    long endTime;
    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

    // -- Construct and store LSH tables --

    construct_lsh_tables(// input
                         vectorLength, numTables, nGroups,
                         nPoints1, points1,
                         nPlanes, hyperplanes,
                         indexGroupMap, indexGroupMapTableLen,
                         // output
                         groupArray,
                         groupSizeMap, groupIndexMap);

#ifdef TIME_LSH
        gettimeofday(&time, NULL);
    endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    printf("Constructing lsh tables:\n");
    printf("   - time: %.3f seconds\n", ((double)endTime - startTime) / 1000);

    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

    // calculate indices into lsh tables for the matching set
    calculate_indexGroupMap(vectorLength, numTables, nPoints2, points2,
                            nPlanes, hyperplanes, indexGroupMapTableLen, indexGroupMap);

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    printf("Calculating group mappings for the other point dataset:\n");
    printf("   - time: %.3f seconds\n", ((double)endTime - startTime) / 1000);

    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

    // Array allocation
    int potentialMatchesMaxLen = nPoints2 * 64;
    int* potentialMatches = (int*)malloc(potentialMatchesMaxLen * sizeof(int));
    int* potentialMatchesIndices = (int*)malloc(nPoints2 * sizeof(int));
    int* potentialMatchesLengths = (int*)malloc(nPoints2 * sizeof(int));

    find_potential_matches(//inputs
                           numTables, nPoints1, nPoints2,
                           indexGroupMap, indexGroupMapTableLen,
                           groupArray, groupSizeMap, groupIndexMap,
                           groupMapTableLen, potentialMatchesMaxLen,
                           // outputs
                           &potentialMatches, potentialMatchesIndices,
                           potentialMatchesLengths);

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    printf("Finding potential matches\n");
    printf("   - time: %.3f seconds\n", ((double)endTime - startTime)/1000);

    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

    // Array allocation and initialization
    // holds the actual matches
    int* lshMatches = (int*)malloc(nPoints2 * sizeof(int));
    memset(lshMatches, -1, nPoints2 * sizeof(int));
    // holds the best match for each point we're finding a match for
    float* bestMatchDists = (float*)malloc(nPoints2 * sizeof(float));
    for (int i = 0; i < nPoints2; i++) {
        bestMatchDists[i] = 1e10;
    }

    // Holds the second best
    int* lshMatches2 = (int*)malloc(nPoints2 * sizeof(int));
    memset(lshMatches2, -1, nPoints2 * sizeof(int));
    // holds the best match for each point we're finding a match for
    float* bestMatchDists2 = (float*)malloc(nPoints2 * sizeof(float));
    for (int i = 0; i < nPoints2; i++) {
        bestMatchDists[i] = 1e10;
    }

    match_points(vectorLength, nPoints2, points2, points1,
                 potentialMatches, potentialMatchesIndices, potentialMatchesLengths,
                 lshMatches, bestMatchDists,
                 lshMatches2, bestMatchDists2);

#ifdef TIME_LSH
    gettimeofday(&time, NULL);
    endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    printf("Matching potential matches\n");
    printf("   - time: %.3f seconds\n", ((double)endTime - startTime) / 1000);
#endif

    // - Check how many matches were correct -
    int correct = 0;
    for (int i = 0; i < nPoints2; i++) {
        correct += lshMatches[i] == i;
    }
    double correctRatio = ((double) correct) / nPoints2;
    printf("Correct ratio: %f\n", correctRatio);

    // // average distances from second best
    // printf("\n");
    // double avgDist1 = 0;
    // double avgDist2 = 0;
    // double avgDist3 = 0;
    // for (int i = 0; i < nPoints2; i++) {
    //     if (bestMatchDists[i] != 1e10 && bestMatchDists2[i] != 1e10) {
    //         double diff = bestMatchDists[i] - bestMatchDists2[i];
    //         diff = diff * diff;
    //         avgDist1 += diff;
    //         avgDist2 += diff * (lshMatches[i] == i);
    //         avgDist3 += diff * (lshMatches[i] != i);
    //     }
    // }
    // avgDist1 /= nPoints2;
    // avgDist2 /= correct;
    // avgDist3 /= (nPoints2 - correct);

    // printf("Average squared distance from second best: %e\n", avgDist1);
    // printf("For correct matches:                       %e\n", avgDist2);
    // printf("For incorrect matches:                     %e\n", avgDist3);
    // printf("(Not counting when only one match was found)\n");

    free(points1);
    free(points2);
    free(hyperplanes);

    free(indexGroupMap);
    free(groupSizeMap);
    free(groupIndexMap);
    free(groupArray);

    free(lshMatches);
    free(bestMatchDists);
    free(lshMatches2);
    free(bestMatchDists2);

    free(potentialMatches);
    free(potentialMatchesIndices);
    free(potentialMatchesLengths);
}

void fill_point_arrays(int nPoints1, int nPoints2, int vectorLength, float noiseScale,
                       float* __restrict__ points1, float* __restrict__ points2)
{
    for (int i = 0; i < ((nPoints1 < nPoints2) ? nPoints1 : nPoints2); i++) {
        float size1 = 0;
        float size2 = 0;
        for (int j = 0; j < vectorLength; j++) {
            // points1[i][:] = random array with values between -1 and 1
            float rn = rand_minus1_to_1();
            points1[i*vectorLength + j] = rn;
            // points2[i][:] = points1[i][:] + noise
            float noise = noiseScale * rand_minus1_to_1();
            points2[i*vectorLength + j] = rn + noise;
            // Keep track of euclidean lengths to normalize after loop
            size1 += points1[i*vectorLength + j] * points1[i*vectorLength + j];
            size2 += points2[i*vectorLength + j] * points2[i*vectorLength + j];
        }

        size1 = sqrt(size1);
        size2 = sqrt(size2);

        // normalize the vectors
        for (int j = 0; j < vectorLength; j++) {
            points1[i*vectorLength + j] /= size1;
            points2[i*vectorLength + j] /= size2;
        }
    }
}

void fill_hyperplanes(int nPlanes, int vectorLength, float* hyperplanes) {

    for (int i = 0; i < nPlanes; i++) {
        float size = 0;
        for (int j = 0; j < vectorLength; j++) {
            // hyperplanes[i][:] = random array with values between -1 and 1
            hyperplanes[i*vectorLength + j] = rand_minus1_to_1();
            // Keep track of euclidean lengths to normalize after loop
            size += hyperplanes[i * vectorLength + j] * hyperplanes[i * vectorLength + j];
        }

        
        size = sqrt(size);
        float invsize = 1/size;
        // normalize the vectors
        for (int j = 0; j < vectorLength; j++) {
            hyperplanes[i * vectorLength + j] *= invsize;
        }
    }
}

void calculate_hash_values(
    // input
    int vectorLength,
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
    // preallocated temporary storage
    int* __restrict__ groupIndexMapTails,
    // output
    int* __restrict__ groupArray,
    int* __restrict__ groupSizeMap,
    int* __restrict__ groupIndexMap)
{
    // find the size of each group
    memset(groupSizeMap, 0, nGroups * sizeof(int));
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
    }

    // Fill groupArray
    for (int i = 0; i < nPoints; i++) {
        // what group did the point get hashed into
        int hashcode = indexGroupMap[i];
        // start index of that group in groupArray + number of elements currently inserted
        int idx = groupIndexMapTails[hashcode]++;
        // add the point to the group
        groupArray[idx] = i;
    }
}

void calculate_indexGroupMap(
    // input
    int vectorLength, int numTables,
    int nPoints, float* __restrict__ points,
    int nPlanes, float* __restrict__ hyperplanes,
    int indexGroupMapTableLen,
    // output
    int* __restrict__ indexGroupMap)
{
    const int hyperplanesTableLen = nPlanes * vectorLength;
    for (int table = 0; table < numTables; table++) {
        float* hyperplanes2 = hyperplanes + table * hyperplanesTableLen;
        int* indexGroupMap2 = indexGroupMap + table * indexGroupMapTableLen;
        // - Match points -
        calculate_hash_values(vectorLength, nPoints, points,
                                nPlanes, hyperplanes2, indexGroupMap2);
    }
}


void construct_lsh_tables(  // input
    int vectorLength, int numTables, int nGroups,
    int nPoints, float* __restrict__ points,
    int nPlanes, float* __restrict__ hyperplanes,
    int* __restrict__ indexGroupMap, int indexGroupMapTableLen,
    // output
    int* __restrict__ groupArray,
    int* __restrict__ groupSizeMap,
    int* __restrict__ groupIndexMap)
{
    const int groupMapTableLen = nGroups;
    const int groupMapLen = numTables * groupMapTableLen;
    const int groupArrayTableLen = nPoints;

    // calculate indices into lsh tables for the original set
    calculate_indexGroupMap(vectorLength, numTables, nPoints, points,
                           nPlanes, hyperplanes, indexGroupMapTableLen, indexGroupMap);

    // temporary values when creating groupIndexMap
    int* groupIndexMapTails = (int*)malloc(groupMapLen * sizeof(int));

    for (int table = 0; table < numTables; table++) {
        int* indexGroupMap2 = indexGroupMap + table * indexGroupMapTableLen;
        int* groupSizeMap2 = groupSizeMap + table * groupMapTableLen;
        int* groupIndexMap2 = groupIndexMap + table * groupMapTableLen;
        int* groupIndexMapTails2 = groupIndexMapTails + table * groupMapTableLen;
        int* groupArray2 = groupArray + table * groupArrayTableLen;

        // - Organize points so they can be indexed by their hash values -
        organize_points_into_groups(nPoints, nGroups, indexGroupMap2,
                                    groupIndexMapTails2,
                                    groupArray2, groupSizeMap2, groupIndexMap2);
    }

    free(groupIndexMapTails);
}

int find_potential_matches(
    // inputs
    int numTables, int nPoints1, int nPoints2,
    int* __restrict__ indexGroupMap, int indexGroupMapTableLen,
    int* __restrict__ groupArray, int* __restrict__ groupSizeMap,
    int* __restrict__ groupIndexMap,
    int groupMapTableLen, int potentialMatchesMaxLen,
    // outputs
    int** __restrict__ potentialMatches,
    int* __restrict__ potentialMatchesIndices,
    int* __restrict__ potentialMatchesLengths)
{
    const int groupArrayTableLen = nPoints1;

    int* checkedArr = (int*)malloc(nPoints1 * sizeof(int));
    memset(checkedArr, -1, nPoints1 * sizeof(int));
    int totalMatchCount = 0;

    // find possible matches
    for (int i = 0; i < nPoints2; i++) {
        potentialMatchesIndices[i] = totalMatchCount;

        for (int table = 0; table < numTables; table++) {
            int* indexGroupMap2 = indexGroupMap + table * indexGroupMapTableLen;
            int* groupSizeMap2 = groupSizeMap + table * groupMapTableLen;
            int* groupIndexMap2 = groupIndexMap + table * groupMapTableLen;
            int* groupArray2 = groupArray + table * groupArrayTableLen;

            // Find the group of elements from groupArray to match with
            int hashcode = indexGroupMap2[i];
            int size = groupSizeMap2[hashcode];
            int startIdx = groupIndexMap2[hashcode];

            // Add points to potentialMatches
            for (int j = startIdx; j < startIdx + size; j++) {
                int idx = groupArray2[j];
                if (checkedArr[idx] != i) {
                    checkedArr[idx] = i;

                    if (totalMatchCount >= potentialMatchesMaxLen) {
                        potentialMatchesMaxLen =
                            double_int_arr_size(potentialMatches, potentialMatchesMaxLen);
                        // printf(" maxLen = %d\n", potentialMatchesMaxLen);
                    }

                    (*potentialMatches)[totalMatchCount] = idx;
                    totalMatchCount++;
                }
            }
        }

        potentialMatchesLengths[i] = totalMatchCount - potentialMatchesIndices[i];
    }

    // printf("\nmatches found: %d\n\n", totalMatchCount);

    free(checkedArr);

    return totalMatchCount;
}

void match_points(
    // inputs
    int vectorLength, int nPoints2,
    float* __restrict__ points2, float* __restrict__ points1,
    int* __restrict__ potentialMatches, int* __restrict__ potentialMatchesIndices,
    int* __restrict__ potentialMatchesLengths,
    // outputs
    int* __restrict__ lshMatches, float* __restrict__ bestMatchDists,
    int* __restrict__ lshMatches2, float* __restrict__ bestMatchDists2)
{
    // printf("new way:");
    for (int i = 0; i < nPoints2; i++) {
        bool changed = false;
        float bestMatchDist = bestMatchDists[i];
        float bestMatchDist2 = bestMatchDist;
        int match = -1;
        int match2 = -1;
        int jStart = potentialMatchesIndices[i];
        int jEnd = potentialMatchesIndices[i] + potentialMatchesLengths[i];

        // printf("i=%d\nidx = \n", i);
        // Match the points
        for (int j = jStart; j < jEnd; j++) {
            int idx = potentialMatches[j];
            // printf(" - %d\n", idx);
            // diff = sum((points2[i][:] - points1[idx][:]) .^ 2);
            float diff = 0;
            for (int k = 0; k < vectorLength; k++) {
                float tmp = points2[i * vectorLength + k] - points1[idx * vectorLength + k];
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

        if (changed) {
            lshMatches[i] = match;
            bestMatchDists[i] = bestMatchDist;
            lshMatches2[i] = match2;
            bestMatchDists2[i] = bestMatchDist2;
        }
    }
}
