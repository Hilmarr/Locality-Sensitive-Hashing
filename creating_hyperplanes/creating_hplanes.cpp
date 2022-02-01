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
    // int nVecs = 200;

    int vectorLength = 128;

    int maxPlanes = 32;
    double* hyperplanes = (double*) malloc(maxPlanes * vectorLength * sizeof(double));

    // number of alternative vectors to choose from in each iteration
    int nAlts = 10;
    double* alts = (double*) malloc(2 * nAlts * vectorLength * sizeof(double));

    double* altsSums = (double*) malloc(2 * nAlts * sizeof(double));
    int* positiveSideCnt = (int*) malloc(nAlts * sizeof(int));
    double* positiveSide = (double*)malloc(nAlts * sizeof(double));
    double* squaredDists = (double*) malloc(nAlts * sizeof(double));

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

            // In the unlikely event where the sum is almost zero, repeat iteration
            if (i > nAlts && sum < 1e-10 && sum > -1e-10) {
                fprintf(stderr,"Unlikely scenario, have to repeat an iteration\n");
                i--;
            }
            // In the even more unlikely event where the size is almost zero,
            // repeat iteration
            if (vecSize < 1e-10 && vecSize > -1e-10) {
                fprintf(stderr, "Unlikely scenario, have to repeat an iteration\n");
                i--;
            }
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
    free(positiveSideCnt);
    free(positiveSide);
    free(squaredDists);

    free(hyperplanes);
}