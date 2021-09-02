#include <iostream>
#include <cmath>

const int nPoints = 1000;
const int vectorLength = 128;
const float noiseScale = 0.1;

// points to be matched
double points1[nPoints][vectorLength];
double points2[nPoints][vectorLength];

inline double rand_minus1_to_1() {
    return ((double)rand() / (RAND_MAX/2)) - 1;
}

void fill_point_arrays();

int main() {
    
    fill_point_arrays();

    double size = 0;
    for (int j = 0; j < vectorLength; j++) {
        size += points1[0][j] * points1[0][j];
        printf("%f\n", points1[0][j]);
    }
    printf("\nsize: %f\n", sqrt(size));
}

void fill_point_arrays() {
    for (int i = 0; i < nPoints; i++) {
        double size1 = 0;
        double size2 = 0;
        for (int j = 0; j < vectorLength; j++) {
            // points1[i] = random array with values between -1 and 1
            double rn = rand_minus1_to_1();
            points1[i][j] = rn;
            // points2[i] = points1[i] + noise
            double noise = noiseScale * rand_minus1_to_1();
            points2[i][j] = rn + noise;
            // Keep track of euclidean lengths to normalize after loop
            size1 += points1[i][j] * points1[i][j];
            size2 += points2[i][j] * points2[i][j];
        }

        size1 = sqrt(size1);
        size2 = sqrt(size2);

        // normalize the vectors
        for (int j = 0; j < vectorLength; j++) {
            points1[i][j] /= size1;
            points2[i][j] /= size2;
        }
    }
}