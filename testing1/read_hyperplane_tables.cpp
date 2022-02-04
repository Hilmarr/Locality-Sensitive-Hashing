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
    


    fclose(fp);
    free(hyperplanes);
}