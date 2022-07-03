#include <iostream>

int main() {
    // char fPath[] = "../hyperplanes.dat";
    // char fPath[] = "../creating_hyperplanes/hyperplanes_attempt1.dat";
    char fPath[] = "../creating_hyperplanes/hyperplanes.dat";
    FILE* fp = fopen(fPath, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error opening %s for reading\n", fPath);
        return -1;
    }

    const int nHplanes = 32;
    const int hPlaneLen = 128;

    float* hPlanes = (float*) malloc(nHplanes * hPlaneLen * sizeof(float));
    int ret = fread(hPlanes, 4, nHplanes*hPlaneLen, fp);
    if (ret != nHplanes*hPlaneLen) {
        printf("Only %d of %d elements read\n", ret, nHplanes*hPlaneLen);
    }

    // int ret = 1;
    // while (ret != 0) {
    //     // move component into memory
    //     float tmp = 0.0;
    //     ret = fread(&tmp, 4, 1, fp);
    //     if (ret != 0) printf("%f ", tmp);
    // }
    // printf("\n");

    // for (int j = 0; j < hPlaneLen; j++) {
    //     printf("%f\n", hPlanes[j]);
    // }

    for (int i = 0; i < nHplanes; i++) {
        for (int j = 0; j < hPlaneLen; j++) {
            printf("%f ", hPlanes[i*hPlaneLen + j]);
        }
        printf("\n\n");
    }

    int cnt = 0;
    for (int i = 0; i < nHplanes; i++) {
        for (int j = i+1; j < nHplanes; j++) {
            double mul = 0;
            for (int k = 0; k < hPlaneLen; k++) {
                mul += hPlanes[i*hPlaneLen + k] * hPlanes[j*hPlaneLen + k];
            }
            // printf("i=%d, j=%d, %f\n", i, j, mul);
            if (mul < 0) mul = -mul;
            if (mul > 1e-6) cnt++;
        }
    }
    printf("Number of non-orthogonal hyperplane pairs = %d\n", cnt);

    cnt = 0;
    for (int i = 0; i < nHplanes; i++) {
        double sum = 0;
        for (int j = 0; j < hPlaneLen; j++) {
            sum += hPlanes[i*hPlaneLen + j];
        }
        // printf("sum(hyperplane[%d]) = %f\n", i, sum);
        if (sum < 0) sum = -sum;
        if (sum > 1e-6) cnt++;
    }
    printf("Number of sums that don't sum to zero = %d\n", cnt);





    fclose(fp);
    free(hPlanes);
}

