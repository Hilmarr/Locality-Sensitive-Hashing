#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>

inline double rand_minus1_to_1() {
    return ((double)rand() / (RAND_MAX/2)) - 1;
}

void fill_point_arrays(int nPoints, int vectorLength, double noiseScale,
                       double* points1, double* points2);

void fill_hyperplanes(int nPlanes, int vectorLength, double* hyperplanes);

int main() {
    const int nPoints = 1000;     // number of points in the first dataset
    const int nPoints2 = nPoints; // number of points in the second dataset
    const int vectorLength = 128;
    const double noiseScale = 0.1;

    // points to be matched
    // LATER: because of their massive size arrays should be dynamically
    //        allocated using 'malloc'
    double points1[nPoints * vectorLength];
    double points2[nPoints2 * vectorLength];

    const int nPlanes = (int)log2(nPoints);

    fill_point_arrays(nPoints, vectorLength, noiseScale, points1, points2);

    // allocate hyperplanes
    double hyperplanes[nPlanes * vectorLength];
    // create normalized hyperplanes represented by vectors of euclidean size 1
    fill_hyperplanes(nPlanes, vectorLength, hyperplanes);

    //  - Arrays to organize groups -
    // LATER: The arrays that whose elements need to
    //        be set to zero are only 'groupSizeMap'
    const int nBoxes = 1 << nPlanes;
    int indexGroupMap[nPoints];  // The group that each point falls into
    int groupSizeMap[nBoxes];    // Amount of points that fall into each group
    int groupIndexMap[nBoxes];   // The starting index for each group in groupArray
    int groupIndexMapTails[nBoxes]; // temporary values

    int groupArray[nPoints]; // Actual group
    int lshMatches[nPoints2]; // holds the actual matches

    // - Calculate hash values, keep track of sizes of each group -
    for (int i = 0; i < nPoints; i++) {
        double* point = &points1[i * vectorLength];
        int hashcode = 0;  //  hashcode will be the group index
        // calculate hash value of the i'th point, store resut in indexGroupMap
        double* hplane = hyperplanes; // first hyperplane
        for (int j = 0; j < nPlanes; j++) {
            // calculate point * hplane
            double vecMul = 0;
            for (int k = 0; k < vectorLength; k++) {
                vecMul += point[k]*hplane[k];
            }
            // set i'th bit to one if point is on "positive" side of hyperplane
            if  (vecMul > 0) {
                hashcode = hashcode | (1 << j);
            }
            // next hyperplane
            hplane += vectorLength;
        }
        indexGroupMap[i] = hashcode;  // save the hashcode
        groupSizeMap[hashcode]++;     // increment the size of the group
    }

    // - Organize points so they can be indexed by their hash values -

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

    // - Match points -

    // first calculate hash values (overwriting values in indexGroupMap)
    for (int i = 0; i < nPoints2; i++) {
        double* point = &points2[i * vectorLength];
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

    // then search
    for (int i = 0; i < nPoints2; i++) {
        // Find the group of elements from groupArray to match with
        double bestFitDist = 1e10;
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
                double tmp = points2[i*vectorLength + k] - points1[idx*vectorLength + k];
                diff += tmp*tmp;
            }
            // check if distance is the lowest distance so far
            if (diff < bestFitDist) {
                // printf("bestFitDist=%.3f   diff=%.3f   match=%d", bestFitDist, diff, match);
                bestFitDist = diff;
                match = idx;
            }
        }

        lshMatches[i] = match;
    }


    // Check how many matches were correct
    int correct = 0;
    for (int i = 0; i < nPoints2; i++) {
        correct += lshMatches[i] == i;
    }
    double  correctRatio = ((double) correct) / nPoints2;
    printf("Correct ratio: %f\n", correctRatio);

    // for (int i = 0; i < nPoints2; i++) {
    //     printf("%d\n", lshMatches[i]);
    // }
}

void fill_point_arrays(int nPoints, int vectorLength, double noiseScale,
                       double* points1, double* points2)
{
    for (int i = 0; i < nPoints; i++) {
        double size1 = 0;
        double size2 = 0;
        for (int j = 0; j < vectorLength; j++) {
            // points1[i][:] = random array with values between -1 and 1
            double rn = rand_minus1_to_1();
            points1[i*vectorLength + j] = rn;
            // points2[i][:] = points1[i][:] + noise
            double noise = noiseScale * rand_minus1_to_1();
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

void fill_hyperplanes(int nPlanes, int vectorLength, double* hyperplanes) {

    for (int i = 0; i < nPlanes; i++) {
        double size = 0;
        for (int j = 0; j < vectorLength; j++) {
            // hyperplanes[i][:] = random array with values between -1 and 1
            hyperplanes[i*vectorLength + j] = rand_minus1_to_1();
            // Keep track of euclidean lengths to normalize after loop
            size += hyperplanes[i * vectorLength + j] * hyperplanes[i * vectorLength + j];
        }

        
        size = sqrt(size);
        double invsize = 1/size;
        // normalize the vectors
        for (int j = 0; j < vectorLength; j++) {
            hyperplanes[i * vectorLength + j] *= invsize;
        }
    }
}