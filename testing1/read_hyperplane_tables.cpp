#include <iostream>

int main() {
    int nTables = 30;
    int nPlanes = 20;
    int vectorLength = 128;

    char fPath[] = "../creating_hyperplanes/hyperplaneTables_100x32x128.dat";
    FILE* fp = fopen(fPath, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error opening %s for reading\n", fPath);
        return -1;
    }

    float* hyperplanes = (float*)malloc(nTables * nPlanes * vectorLength * sizeof(float));

    for (int t = 0; t < nTables; t++) {
        int ret = fread(&hyperplanes[t*nPlanes*vectorLength], 4, nPlanes * vectorLength, fp);
        if (ret != nPlanes * vectorLength) {
            fprintf(stderr, "Table: %d, Only %d of %d elements read\n", t, ret, nPlanes * vectorLength);
            break;
        }
        fseek(fp, (32-nPlanes)*128*4, SEEK_CUR);
    }

    for (int t = 0; t < nTables; t++) {
        float* hp = &hyperplanes[t*nPlanes*vectorLength];

        int cnt = 0;
        for (int i = 0; i < nPlanes; i++) {
            for (int j = i + 1; j < nPlanes; j++) {
                double mul = 0;
                for (int k = 0; k < vectorLength; k++) {
                    mul += hp[i * vectorLength + k] * hp[j * vectorLength + k];
                }
                // printf("i=%d, j=%d, %f\n", i, j, mul);
                if (mul < 0) mul = -mul;
                if (mul > 1e-6) cnt++;
            }
        }
        printf("Number of non-orthogonal hyperplane pairs = %d\n", cnt);

        cnt = 0;
        for (int i = 0; i < nPlanes; i++) {
            double sum = 0;
            for (int j = 0; j < vectorLength; j++) {
                sum += hp[i * vectorLength + j];
            }
            // printf("sum(hyperplane[%d]) = %f\n", i, sum);
            if (sum < 0) sum = -sum;
            if (sum > 1e-6) cnt++;
        }
        printf("Number of sums that don't sum to zero = %d\n", cnt);
    }

    fclose(fp);
    free(hyperplanes);
}