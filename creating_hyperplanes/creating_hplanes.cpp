#include <cmath>
#include <iostream>

#include "read_vector_file_funcs.h"

inline double rand_minus1_to_1() {
    return ((double)rand() / (RAND_MAX / 2)) - 1;
}

int main(int argc, char **argv) {
    srand(0); // Stable results for now, comment out later

    // Check number of arguments
    if (argc < 2) {
        fprintf(stderr, "Not enough arguments\n");
        fprintf(stderr, "You need to also provide a file with learning vectors\n");
        return -1;
    }

    int* tmp1;
    int nVecsRead = read_vector_file2(argv[1], &tmp1) / 128;
    float* vecs = (float*)tmp1;
    // number of vectors we actually use
    int nVecs = 10;

    int vectorLength = 128;

    int maxPlanes = 3;
    double* hyperplanes = (double*) malloc(maxPlanes * vectorLength * sizeof(double));

    // number of alternative vectors to choose from in each iteration
    int nAlts = 10;
    double* alts = (double*) malloc(2 * nAlts * vectorLength * sizeof(double));

    double* altsSums = (double*) malloc(2 * nAlts * sizeof(double));
    double* positiveSide = (double*)malloc(nAlts * sizeof(double));
    double* squaredDists = (double*) malloc(nAlts * sizeof(double));
    // double* costs = (double*)malloc(nAlts * sizeof(double));

    for (int nPlanes = 0; nPlanes < maxPlanes; nPlanes++) {
        // Create possible hyperperplanes + some extra used for adjustments
        // also calculate sums
        for (int i = 0; i < 2 * nAlts; i++) {
            // Initialize vector as random vector
            for (int j = 0; j < vectorLength; j++) {
                alts[i*vectorLength + j] = rand_minus1_to_1();
            }
            
            // Make vector orthogonal to hyperplanes
            for (int k = 0; k < nPlanes; k++) {
                // Calculate dot product
                double prod = 0;
                for (int j = 0; j < vectorLength; j++) {
                    prod += alts[i*vectorLength + j] * hyperplanes[k*vectorLength + j];
                }
                // Remove projection of vector onto hyperplane from the vector
                for (int j = 0; j < vectorLength; j++) {
                    alts[i*vectorLength + j] -= prod * hyperplanes[k*vectorLength + j];
                }
            }

            // Calculate sum
            double sum = 0;
            double vecSize = 0;
            for (int j = 0; j < vectorLength; j++) {
                double val = alts[i*vectorLength + j]; 
                sum += val;
                vecSize += val*val;
            }
            altsSums[i] = sum;

            // In the unlikely event where the sum is almost zero: repeat iteration
            if (i > nAlts && sum < 1e-10 && sum > -1e-10) {
                fprintf(stderr,"Unlikely scenario occurred, have to repeat an iteration\n");
                i--;
                continue;
            }
            // In the unlikely event where the squared euclidean length
            // is almost zero: repeat iteration
            if (vecSize < 1e-10 && vecSize > -1e-10) {
                fprintf(stderr, "Unlikely scenario occurred, have to repeat an iteration\n");
                i--;
                continue;
            }
        }

        // Make alternative vectors sum to zero
        for (int i = 0; i < nAlts; i++) {
            double factor = altsSums[i] / altsSums[nAlts+i];
            printf("altsSums[%d]=%.3f, altsSums[%d]=%.3f, factor=%f\n",
                   i, altsSums[i], nAlts + i, altsSums[nAlts + i], factor);
            for (int j = 0; j < vectorLength; j++) {
                alts[i*vectorLength + j] -= factor
                    * alts[(nAlts+i)*vectorLength + j];
            }
        }

        bool repeat = false;

        // Normalize vectors
        for (int i = 0; i < nAlts; i++) {
            // Find euclidean length
            double vecSize = 0;
            for (int j = 0; j < vectorLength; j++) {
                double val = alts[i * vectorLength + j];
                vecSize += val * val;
            }
            vecSize = sqrt(vecSize);
            // Check that each hyperplane still has a euclidean length bigger than
            // zero, if not, repeat outer loop
            if (vecSize < 1e-10 && vecSize > -1e-10) {
                repeat = true;
                break;
            }
            // Normalize
            for (int j = 0; j < vectorLength; j++) {
                alts[i * vectorLength + j] *= vecSize;
            }
        }

        // Since error will never happen in practice, all the program will do is just
        // redo the iteration.
        if (repeat) {
            fprintf(stderr, "Unlikely scenario occurred, have to repeat outer loop iteration\n");
            nPlanes--;
            continue;
        }

        // Find ratio that falls on the positive vs negative side,
        // and distances from hyperplanes
        for (int i = 0; i < nAlts; i++) {
            int positiveCnt = 0;
            double squaredDistsAcc = 0;
            for (int k = 0; k < nVecs; k++) {
                double prod = 0;
                for (int j = 0; j < vectorLength; j++) {
                    prod += alts[i*vectorLength+j] + vecs[k*vectorLength + j];
                }
                squaredDistsAcc += prod*prod;
                positiveCnt += prod > 0;
            }
            
            double tmp = (double) positiveCnt;
            tmp = 0.5 - (tmp / nVecs);
            tmp = tmp*tmp;
            positiveSide[i] = tmp;

            squaredDists[i] = 1 / (1 + (squaredDistsAcc / nVecs));
        }

        // Calculate costs for each potential hyperplane
        double acc1 = 0;
        double acc2 = 0;
        for (int i = 0; i < nAlts; i++) {
            acc1 = positiveSide[i] * positiveSide[i];
            acc2 = squaredDists[i] * squaredDists[i];
        }
        acc1 /= nAlts;
        acc2 /= nAlts;
        for (int i = 0; i < nAlts; i++) {
            positiveSide[i] /= acc1;
            squaredDists[i] /= acc2;
        }
        double minCost = positiveSide[0] + squaredDists[0];
        int minCostIdx = 0;
        for (int i = 1; i < nAlts; i++) {
            double cost = positiveSide[i] + squaredDists[i];
            if (cost < minCost) {
                minCost = cost;
                minCostIdx = i;
            }
        }

        // Choose hyperplane with the lowest cost
        for (int i = 0; i < vectorLength; i++) {
            hyperplanes[nPlanes*vectorLength + i] =  alts[minCostIdx*vectorLength + i];
        }


        // // Calculate sums
        // for (int i = 0; i < 2 * nAlts; i++) {
        //     double sum = 0;
        //     for (int j = 0; j < vectorLength; j++) {
        //         sum += alts[i*vectorLength + j];
        //     }
        //     // The sum of 128 random numbers should never be zero, but if it is
        //     // we exit the program
        //     // if (sum < 1e-12) {
        //     //     fprintf(stderr, "Very unlikely critical error, exiting program\n");
        //     //     return -1;
        //     // }
        //     printf("sum = %f\n", sum);
        //     altsSums[i] = sum;
        // }

        // // Make alternative vectors to choose from orthogonal
        // for (int i = 0; i < nAlts; i++) {
        //     double factor = altsSums[i] / altsSums[i+nAlts];
            
        // }
    }

    free(vecs);

    free(alts);

    free(altsSums);
    free(positiveSide);
    free(squaredDists);
    // free(costs);

    free(hyperplanes);
}