#include "read_vector_file_funcs.h"
#include <iostream>

//  -- Taking care to not allocate more memory than necessary --
int read_vector_file(char* fPath, int** arr) {
    FILE* fp = fopen(fPath, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error opening %s for reading\n", fPath);
        return -1;
    }

    // find the size of the file in terms of elements
    int nElements = 0;
    int ret = 1;
    while (ret != 0) {
        int d = 0;
        ret = fread(&d, 4, 1, fp);

        nElements += d;
        fseek(fp, d * 4, SEEK_CUR);

        // printf("%d\n", d);
    }
    rewind(fp);

    // printf("nElements: %d\n", nElements);

    // where the data should be loaded into
    int* data = (int*)malloc(nElements * 4);

    // copy the data from the file into memory
    int* data_ptr = data;
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

    *arr = data;

    fclose(fp);

    return nElements;
}

// -- Allocating enough memory for the entire file --
int read_vector_file2(char* fPath, int** arr) {
    FILE* fp = fopen(fPath, "rb");
    if (fp == NULL) {
        fprintf(stderr, "Error opening %s for reading\n", fPath);
        return -1;
    }

    fseek(fp, 0L, SEEK_END);
    int filesize = ftell(fp);
    rewind(fp);

    // where the data should be loaded into
    int* data = (int*)malloc(filesize);

    int nElements = 0;

    // copy the data from the file into memory
    int* data_ptr = data;
    int ret = 1;
    while (ret != 0) {
        // read nr of elements in component
        int d = 0;
        ret = fread(&d, 4, 1, fp);
        nElements += d;
        if (ret == 0) break;
        // move component into memory
        ret = fread(data_ptr, 4, d, fp);
        data_ptr += d;
    }

    *arr = data;

    fclose(fp);

    return nElements;
}