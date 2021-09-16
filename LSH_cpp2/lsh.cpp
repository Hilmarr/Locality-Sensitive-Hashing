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
 * @param nPoints1 : Length of points1 array
 * @param nPoints2 : Length of points2 array
 * @param vectorLength : Number of dimensions for each point
 * @param noisScale : The amount of noise relative to the length of each dimension
 * @param points1 : Output, point array 1
 * @param points2 : Output, point array 2
 */
void fill_point_arrays(int nPoints1, int nPoints2, int vectorLength, float noiseScale,
                       float* points1, float* points2);

/**
 * Fills the hyperplane array with random hyperplane vectors.
 * @param nPlanes : Number of hyperplanes
 * @param vectorLength : Number of dimensions for each hyperplane
 * @param hyperplanes: Output, the returned hyperplanes
 */
void fill_hyperplanes(int nPlanes, int vectorLength, float* hyperplanes);

void calculate_hash_values(// input
                           int vectorLength,
                           int nPoints, float* points,
                           int nPlanes, float* hyperplanes,
                           // output
                           int* indexGroupMap);

void organize_points_into_groups(// input
                                 int nPoints, int nBoxes,
                                 int* indexGroupMap,
                                 // preallocated temporary storage
                                 int* groupIndexMapTails,
                                 // output
                                 int* groupSizeMap, int* groupIndexMap, int* groupArray);

void calculate_indexGroupMap(// input
                             int vectorLength, int numTables,
                             int nPoints, float* points,
                             int nPlanes, float* hyperplanes,
                             int hyperplanesTableLen, int indexGroupMapTableLen,
                             // output
                             int* indexGroupMap);

void construct_lsh_tables(// input
                          int vectorLength, int numTables, int nBoxes,
                          int nPoints, float* points,
                          int nPlanes, float* hyperplanes, int hyperplanesTableLen,
                          int* indexGroupMap, int indexGroupMapTableLen,
                          // output
                          int* groupArray, int groupArrayTableLen,
                          int* groupSizeMap, int* groupIndexMap);

int find_potential_matches(// inputs
                           int vectorLength, int numTables, int nPoints1, int nPoints2,
                           int* indexGroupMap, int indexGroupMapTableLen,
                           int* groupSizeMap, int* groupIndexMap, int groupMapTableLen, 
                           int* groupArray, int groupArrayTableLen,
                           int potentialMatchesMaxLen,
                           // outputs
                           int* potentialMatches, int* potentialMatchesIndices,
                           int* potentialMatchesLengths);

void match_points(// inputs
                  int vectorLength, int nPoints2, float* points2, float* points1,
                  int* potentialMatches, int* potentialMatchesIndices, int* potentialMatchesLengths,
                  // outputs
                  int* lshMatches, float* bestMatchDists);


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
    const int nPoints1 = 10000;     // number of points in the first dataset
    const int nPoints2 = nPoints1;  // number of points in the second dataset
    const int vectorLength = 128;
    const float noiseScale = 0.3;
    const int numTables = 8;

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

    const int nBoxes = 1 << nPlanes;
    // the group that each point falls into
    const int indexGroupMapTableLen = ((nPoints1 > nPoints2) ? nPoints1 : nPoints2) * sizeof(int);
    const int indexGroupMapLen = numTables * indexGroupMapTableLen;
    int* indexGroupMap = (int*)malloc(indexGroupMapLen * sizeof(int));

    // in case we're changing sizes
    const int groupMapTableLen = nBoxes;
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

    // -- Construct and store LSH tables using

    construct_lsh_tables(// input
                         vectorLength, numTables, nBoxes,
                         nPoints1, points1,
                         nPlanes, hyperplanes, hyperplanesTableLen,
                         indexGroupMap, indexGroupMapTableLen,
                         // output
                         groupArray, groupArrayTableLen,
                         groupSizeMap, groupIndexMap);

    // calculate indices into lsh tables for the matching set
    calculate_indexGroupMap(vectorLength, numTables, nPoints2, points2,
                            nPlanes, hyperplanes, hyperplanesTableLen,
                            indexGroupMapTableLen, indexGroupMap);

#ifdef TIME_LSH
        gettimeofday(&time, NULL);
    endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    printf("Constructing lsh tables:\n");
    printf("   - time: %.3f seconds\n", ((double)endTime - startTime) / 1000);

    gettimeofday(&time, NULL);
    startTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif

    // Array allocation
    int totalMatchCount = 0;
    int potentialMatchesMaxLen = nPoints2 * 64;
    int* potentialMatches = (int*)malloc(potentialMatchesMaxLen * sizeof(int));
    int* potentialMatchesIndices = (int*)malloc(nPoints2 * sizeof(int));
    int* potentialMatchesLengths = (int*)malloc(nPoints2 * sizeof(int));

    find_potential_matches(//inputs
                           vectorLength, numTables, nPoints1, nPoints2,
                           indexGroupMap, indexGroupMapTableLen,
                           groupSizeMap, groupIndexMap, groupMapTableLen, 
                           groupArray, groupArrayTableLen,
                           potentialMatchesMaxLen,
                           // outputs
                           potentialMatches, potentialMatchesIndices,
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

    match_points(vectorLength, nPoints2, points2, points1,
                 potentialMatches, potentialMatchesIndices, potentialMatchesLengths,
                 lshMatches, bestMatchDists);

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
    double  correctRatio = ((double) correct) / nPoints2;
    printf("Correct ratio: %f\n", correctRatio);

    free(points1);
    free(points2);
    free(hyperplanes);

    free(indexGroupMap);
    free(groupSizeMap);
    free(groupIndexMap);
    free(groupArray);
    free(lshMatches);
    free(bestMatchDists);

    free(potentialMatches);
    free(potentialMatchesIndices);
    free(potentialMatchesLengths);
}

void fill_point_arrays(int nPoints1, int nPoints2, int vectorLength, float noiseScale,
                       float* points1, float* points2)
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

void calculate_hash_values(// input
                           int vectorLength,
                           int nPoints, float* points,
                           int nPlanes, float* hyperplanes,
                           // output
                           int* indexGroupMap)
{
    // - Calculate hash values, keep track of sizes of each group -
    for (int i = 0; i < nPoints; i++) {
        float* point = &points[i * vectorLength];
        int hashcode = 0;  //  hashcode will be the group index
        // calculate hash value of the i'th point, store resut in indexGroupMap
        float* hplane = hyperplanes;  // first hyperplane
        for (int j = 0; j < nPlanes; j++) {
            // calculate point * hplane
            float vecMul = 0;
            for (int k = 0; k < vectorLength; k++) {
                vecMul += point[k] * hplane[k];
            }
            // set i'th bit to one if point is on "positive" side of hyperplane
            if (vecMul > 0) {
                hashcode = hashcode | (1 << j);
            }
            // next hyperplane
            hplane += vectorLength;
        }
        indexGroupMap[i] = hashcode;  // save the hashcode
    }
}

void organize_points_into_groups(// input
                                 int nPoints, int nBoxes,
                                 int* indexGroupMap,
                                 // preallocated storage, used to store temproary data
                                 int* groupIndexMapTails,
                                 // output
                                 int* groupSizeMap, int* groupIndexMap, int* groupArray)
{
    // find the size of each group
    memset(groupSizeMap, 0, nBoxes * sizeof(int));
    for (int i = 0; i < nPoints; i++) {
        groupSizeMap[indexGroupMap[i]]++;
    }

    // Find group indices into the group array 'groupArray'
    // // find group indices using exclusive scan of the group sizes
    //  std::exclusive_scan(groupSizeMap, groupSizeMap + nBoxes - 1);
    // doing it manually for now
    int cnt = 0;
    for (int i = 0; i < nBoxes; i++) {
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

void calculate_indexGroupMap(// input
                             int vectorLength, int numTables,
                             int nPoints, float* points,
                             int nPlanes, float* hyperplanes,
                             int hyperplanesTableLen, int indexGroupMapTableLen,
                             // output
                             int* indexGroupMap)
{
    for (int table = 0; table < numTables; table++) {
        float* hyperplanes2 = hyperplanes + table * hyperplanesTableLen;
        int* indexGroupMap2 = indexGroupMap + table * indexGroupMapTableLen;
        // - Match points -
        calculate_hash_values(vectorLength, nPoints, points,
                                nPlanes, hyperplanes2, indexGroupMap2);
    }
}


void construct_lsh_tables(// input
                          int vectorLength, int numTables, int nBoxes,
                          int nPoints, float* points,
                          int nPlanes, float* hyperplanes, int hyperplanesTableLen,
                          int* indexGroupMap, int indexGroupMapTableLen,
                          // output
                          int* groupArray, int groupArrayTableLen,
                          int* groupSizeMap, int* groupIndexMap)
{
    // in case we're changing sizes
    const int groupMapTableLen = nBoxes;
    const int groupMapLen = numTables * groupMapTableLen;

    // calculate indices into lsh tables for the original set
    calculate_indexGroupMap(vectorLength, numTables, nPoints, points,
                           nPlanes, hyperplanes, hyperplanesTableLen,
                           indexGroupMapTableLen, indexGroupMap);

    // temporary values when creating groupIndexMap
    int* groupIndexMapTails = (int*)malloc(groupMapLen * sizeof(int));

    for (int table = 0; table < numTables; table++) {
        int* indexGroupMap2 = indexGroupMap + table * indexGroupMapTableLen;
        int* groupSizeMap2 = groupSizeMap + table * groupMapTableLen;
        int* groupIndexMap2 = groupIndexMap + table * groupMapTableLen;
        int* groupIndexMapTails2 = groupIndexMapTails + table * groupMapTableLen;
        int* groupArray2 = groupArray + table * groupArrayTableLen;

        // - Organize points so they can be indexed by their hash values -
        organize_points_into_groups(nPoints, nBoxes, indexGroupMap2,
                                    groupIndexMapTails2,
                                    groupSizeMap2, groupIndexMap2, groupArray2);
    }

    free(groupIndexMapTails);
}

int find_potential_matches(// inputs
                           int vectorLength, int numTables, int nPoints1, int nPoints2,
                           int* indexGroupMap, int indexGroupMapTableLen,
                           int* groupSizeMap, int* groupIndexMap, int groupMapTableLen, 
                           int* groupArray, int groupArrayTableLen,
                           int potentialMatchesMaxLen,
                           // outputs
                           int* potentialMatches, int* potentialMatchesIndices,
                           int* potentialMatchesLengths)
{
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
                            double_int_arr_size(&potentialMatches, potentialMatchesMaxLen);
                        // printf(" maxLen = %d\n", potentialMatchesMaxLen);
                    }

                    potentialMatches[totalMatchCount] = idx;
                    totalMatchCount++;
                }
            }
        }

        potentialMatchesLengths[i] = totalMatchCount - potentialMatchesIndices[i];
    }

    free(checkedArr);

    return totalMatchCount;
}

void match_points(// inputs
                  int vectorLength, int nPoints2, float* points2, float* points1,
                  int* potentialMatches, int* potentialMatchesIndices, int* potentialMatchesLengths,
                  // outputs
                  int* lshMatches, float* bestMatchDists)
{
    // printf("new way:");
    for (int i = 0; i < nPoints2; i++) {
        bool changed = false;
        float bestFitDist = bestMatchDists[i];
        int match = -1;
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
            if (diff < bestFitDist) {
                // printf("bestFitDist=%.3f   diff=%.3f   match=%d", bestFitDist, diff, match);
                bestFitDist = diff;
                match = idx;
                changed = true;
            }
        }

        if (changed) {
            lshMatches[i] = match;
            bestMatchDists[i] = bestFitDist;
        }
    }
}
