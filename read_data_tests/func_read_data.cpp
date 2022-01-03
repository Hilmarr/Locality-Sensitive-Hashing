#include <iostream>

int main(int argc, char** argv) {
    // Check number of arguments
    if (argc < 2) {
        fprintf(stderr, "Program needs to be given a file as an argument\n");
        return -1;
    }

    // Open file
    char* fPath = argv[1];
    FILE* fp = fopen(fPath, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error opening %s for reading\n", fPath);
        return -1;
    }

    // find the size of the file in terms of elements
    int elements = 0;
    int ret = 1;
    while (ret != 0) {
        int d = 0;
        ret = fread(&d, 4, 1, fp);

        elements += d;
        fseek(fp, d * 4, SEEK_CUR);

        // printf("%d\n", d);
    }
    rewind(fp);

    printf("Elements: %d\n", elements);

    // where the data should be loaded into
    float* data = (float*) malloc(elements * 4);

    // copy the data from the file into memory
    float* data_ptr = data;
    ret = 1;
    while (ret != 0) {
        // read nr of elements in component
        int d = 0;
        ret = fread(&d, 4, 1, fp);
        if (ret == 0) break;
        // move component into memory
        ret = fread(data_ptr, 4, d, fp);
        data_ptr += d;
    }

    for (int i = 0; i < 128; i++) {
        printf("%f\n", data[i]);
    }
}
