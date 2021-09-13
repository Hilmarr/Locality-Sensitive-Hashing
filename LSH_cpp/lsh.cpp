#include <cmath>
#include <cstring>
#include <iostream>
#include <iterator>
#include <numeric>

inline double rand_minus1_to_1() {
    return ((double)rand() / (RAND_MAX / 2)) - 1;
}

void fill_point_arrays(int nPoints1, int nPoints2, int vectorLength, double noiseScale,
                       double* points1, double* points2);

void fill_hyperplanes(int nPlanes, int vectorLength, double* hyperplanes);

void calculate_hash_values(int nPoints, int nPlanes, int vectorLength,
                           double* points, double* hyperplanes,
                           int* indexGroupMap);

void organize_points_into_groups(int nPoints, int nBoxes,
                                 int* indexGroupMap, int* groupSizeMap,
                                 int* groupIndexMap, int* groupIndexMapTails,
                                 int* groupArray);

void lsh_match_points(int nPoints2, int vectorLength, double* points1,
                      double* points2, int* indexGroupMap, int* groupSizeMap,
                      int* groupIndexMap, int* groupArray, int* lshMatches,
                      double* bestMatchDists);

int main() {
    const int nPoints = 10000;     // number of points in the first dataset
    const int nPoints2 = nPoints;  // number of points in the second dataset
    const int vectorLength = 128;
    const double noiseScale = 0.3;
    const int numTables = 8;

    // points to be matched
    double* points1 = (double*)malloc(numTables * nPoints * vectorLength * sizeof(double));
    double* points2 = (double*)malloc(numTables * nPoints2 * vectorLength * sizeof(double));

    const int nPlanes = (int)log2(nPoints);

    fill_point_arrays(nPoints, nPoints2, vectorLength, noiseScale, points1, points2);

    // allocate hyperplanes
    const int hyperplanesTableLen = nPlanes * vectorLength;
    const int hyperplanesLen = numTables * hyperplanesTableLen;
    double* hyperplanes = (double*)malloc(hyperplanesLen * sizeof(double));

    //  - Arrays to organize groups -
    const int nBoxes = 1 << nPlanes;
    // the group that each point falls into
    const int indexGroupMapTableLen = ((nPoints > nPoints2) ? nPoints : nPoints2) * sizeof(int);
    const int indexGroupMapLen = numTables * indexGroupMapTableLen;
    int* indexGroupMap = (int*)malloc(indexGroupMapLen * sizeof(int));

    // in case we're changing sizes
    const int groupMapTableLen = nBoxes;
    const int groupMapLen = numTables * groupMapTableLen;
    // amount of points that fall into each group
    int* groupSizeMap = (int*)calloc(groupMapLen, sizeof(int));
    // the starting index for each group in groupArray
    int* groupIndexMap = (int*)malloc(groupMapLen * sizeof(int));
    // temporary values
    int* groupIndexMapTails = (int*)malloc(groupMapLen * sizeof(int));

    // Actual groups
    const int groupArrayTableLen = nPoints;
    const int groupArrayLen = numTables * groupArrayTableLen;
    int* groupArray = (int*)malloc(groupArrayLen * sizeof(int));

    // holds the actual matches
    int* lshMatches = (int*)malloc(nPoints2 * sizeof(int));
    memset(lshMatches, -1, nPoints2 * sizeof(int));
    // holds the best match for each point we're finding a match for
    double* bestMatchDists = (double*)malloc(nPoints2 * sizeof(double));
    for (int i = 0; i < nPoints2; i++) {
        bestMatchDists[i] = 1e10;
    }

    for (int table = 0; table < numTables; table++) {
        double* hyperplanes2 = hyperplanes + table * hyperplanesTableLen;
        int* indexGroupMap2 = indexGroupMap + table * indexGroupMapTableLen;
        int* groupSizeMap2 = groupSizeMap + table * groupMapTableLen;
        int* groupIndexMap2 = groupIndexMap + table * groupMapTableLen;
        int* groupIndexMapTails2 = groupIndexMapTails + table * groupMapTableLen;
        int* groupArray2 = groupArray + table * groupArrayTableLen;

        // create normalized hyperplanes represented by vectors of euclidean size 1
        fill_hyperplanes(nPlanes, vectorLength, hyperplanes2);

        calculate_hash_values(nPoints, nPlanes, vectorLength,
                              points1, hyperplanes2, indexGroupMap2);

        // - Organize points so they can be indexed by their hash values -
        organize_points_into_groups(nPoints, nBoxes, indexGroupMap2, groupSizeMap2,
                                    groupIndexMap2, groupIndexMapTails2, groupArray2);

        // - Match points -

        calculate_hash_values(nPoints2, nPlanes, vectorLength,
                              points2, hyperplanes2, indexGroupMap2);

        lsh_match_points(nPoints2, vectorLength, points1, points2, indexGroupMap2,
                         groupSizeMap2, groupIndexMap2, groupArray2,
                         lshMatches, bestMatchDists);
    }

    // - Check how many matches were correct -
    int correct = 0;
    for (int i = 0; i < nPoints2; i++) {
        correct += lshMatches[i] == i;
    }
    double correctRatio = ((double)correct) / nPoints2;
    printf("Correct ratio: %f\n", correctRatio);

    free(points1);
    free(points2);
    free(hyperplanes);

    free(indexGroupMap);
    free(groupSizeMap);
    free(groupIndexMap);
    free(groupIndexMapTails);
    free(groupArray);
    free(lshMatches);
    free(bestMatchDists);
}

void fill_point_arrays(int nPoints1, int nPoints2, int vectorLength, double noiseScale,
                       double* points1, double* points2) {
    for (int i = 0; i < ((nPoints1 < nPoints2) ? nPoints1 : nPoints2); i++) {
        double size1 = 0;
        double size2 = 0;
        for (int j = 0; j < vectorLength; j++) {
            // points1[i][:] = random array with values between -1 and 1
            double rn = rand_minus1_to_1();
            points1[i * vectorLength + j] = rn;
            // points2[i][:] = points1[i][:] + noise
            double noise = noiseScale * rand_minus1_to_1();
            points2[i * vectorLength + j] = rn + noise;
            // Keep track of euclidean lengths to normalize after loop
            size1 += points1[i * vectorLength + j] * points1[i * vectorLength + j];
            size2 += points2[i * vectorLength + j] * points2[i * vectorLength + j];
        }

        size1 = sqrt(size1);
        size2 = sqrt(size2);

        // normalize the vectors
        for (int j = 0; j < vectorLength; j++) {
            points1[i * vectorLength + j] /= size1;
            points2[i * vectorLength + j] /= size2;
        }
    }
}

void fill_hyperplanes(int nPlanes, int vectorLength, double* hyperplanes) {
    for (int i = 0; i < nPlanes; i++) {
        double size = 0;
        for (int j = 0; j < vectorLength; j++) {
            // hyperplanes[i][:] = random array with values between -1 and 1
            hyperplanes[i * vectorLength + j] = rand_minus1_to_1();
            // Keep track of euclidean lengths to normalize after loop
            size += hyperplanes[i * vectorLength + j] * hyperplanes[i * vectorLength + j];
        }

        size = sqrt(size);
        double invsize = 1 / size;
        // normalize the vectors
        for (int j = 0; j < vectorLength; j++) {
            hyperplanes[i * vectorLength + j] *= invsize;
        }
    }
}

void calculate_hash_values(int nPoints, int nPlanes, int vectorLength,
                           double* points, double* hyperplanes,
                           int* indexGroupMap) {
    // - Calculate hash values, keep track of sizes of each group -
    for (int i = 0; i < nPoints; i++) {
        double* point = &points[i * vectorLength];
        int hashcode = 0;  //  hashcode will be the group index
        // calculate hash value of the i'th point, store resut in indexGroupMap
        double* hplane = hyperplanes;  // first hyperplane
        for (int j = 0; j < nPlanes; j++) {
            // calculate point * hplane
            double vecMul = 0;
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

void organize_points_into_groups(int nPoints, int nBoxes,
                                 int* indexGroupMap, int* groupSizeMap,
                                 int* groupIndexMap, int* groupIndexMapTails,
                                 int* groupArray) {
    // // find the size of each group
    // memset(groupSizeMap, 0, nBoxes * sizeof(int));
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

void lsh_match_points(int nPoints2, int vectorLength, double* points1,
                      double* points2, int* indexGroupMap, int* groupSizeMap,
                      int* groupIndexMap, int* groupArray, int* lshMatches,
                      double* bestMatchDists) {
    for (int i = 0; i < nPoints2; i++) {
        // Find the group of elements from groupArray to match with
        bool changed = false;
        double bestFitDist = bestMatchDists[i];
        int match = -1;
        int hashcode = indexGroupMap[i];
        int size = groupSizeMap[hashcode];
        int startIdx = groupIndexMap[hashcode];

        // Match the points
        for (int j = startIdx; j < startIdx + size; j++) {
            int idx = groupArray[j];
            // diff = sum((points2[i][:] - points1[idx][:]) .^ 2);
            double diff = 0;
            for (int k = 0; k < vectorLength; k++) {
                double tmp = points2[i * vectorLength + k] - points1[idx * vectorLength + k];
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
