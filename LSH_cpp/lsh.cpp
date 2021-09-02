#include <iostream>
#include <cmath>

inline double rand_minus1_to_1() {
    return ((double)rand() / (RAND_MAX/2)) - 1;
}

void fill_point_arrays(int nPoints, int vectorLength, int noiseScale,
                       double* points1, double* points2);

int main() {
    const int nPoints = 1000;
    const int vectorLength = 128;
    const float noiseScale = 0.1;

    // points to be matched
    double points1[nPoints * vectorLength];
    double points2[nPoints * vectorLength];

    const int nPlanes = (int)log2(nPoints);

    fill_point_arrays(nPoints, vectorLength, noiseScale, points1, points2);

    double size = 0;
    for (int j = 0; j < vectorLength; j++) {
        size += points1[j] * points1[j];
        printf("%f\n", points1[j]);
    }
    printf("\nsize: %f\n", sqrt(size));

    // double hyperplanes[nPlanes][vectorLength];
}

void fill_point_arrays(int nPoints, int vectorLength, int noiseScale,
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