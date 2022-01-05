#include <iostream>

//  -- Taking care to not allocate more memory than necessary --
int read_vector_file(char* fPath, void** arr) {
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
    
    *arr = (void*) data;

    fclose(fp);

    return nElements;
}

// -- Allocating enough memory for the entire file --
int read_vector_file2(char* fPath, void** arr) {
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

    *arr = (void*) data;

    fclose(fp);

    return nElements;
}

int main(int argc, char** argv) {
    // Check number of arguments
    if (argc < 2) {
        fprintf(stderr, "Program needs to be given a file as an argument\n");
        return -1;
    }

    // // Open file
    // char* fPath = argv[1];
    // FILE* fp = fopen(fPath, "rb");
    // if (fp == NULL) {
    //     fprintf(stderr, "Error opening %s for reading\n", fPath);
    //     return -1;
    // }

    // // find the size of the file in terms of elements
    // int nElements = 0;
    // int ret = 1;
    // while (ret != 0) {
    //     int d = 0;
    //     ret = fread(&d, 4, 1, fp);

    //     nElements += d;
    //     fseek(fp, d * 4, SEEK_CUR);

    //     // printf("%d\n", d);
    // }
    // rewind(fp);

    // printf("nElements: %d\n", nElements);

    // // where the data should be loaded into
    // float* data = (float*) malloc(nElements * 4);

    // // copy the data from the file into memory
    // float* data_ptr = data;
    // ret = 1;
    // while (ret != 0) {
    //     // read nr of elements in component
    //     int d = 0;
    //     ret = fread(&d, 4, 1, fp);
    //     if (ret == 0) break;
    //     // move component into memory
    //     ret = fread(data_ptr, 4, d, fp);
    //     data_ptr += d;
    // }

    // for (int i = 0; i < 128; i++) {
    //     printf("%f\n", data[i]);
    // }

    //  -- Taking care to not allocate more memory than necessary --
    
    float* data;
    int nElements;

    nElements = read_vector_file(argv[1], (void**) &data);

    // for (int i = 0; i < 128; i++) {
    //     printf("%f\n", data[i]);
    // }

    printf("1 nElements: %d\n", nElements);
    printf("1 number of vectors: %d\n", nElements / 128);

    free(data);


    // -- Allocating enough memory for the entire file --

    nElements = read_vector_file2(argv[1], (void**)& data);

    // for (int i = 0; i < 128; i++) {
    //     printf("%f\n", data[i]);
    // }

    printf("2 nElements: %d\n", nElements);
    printf("2 number of vectors: %d\n", nElements / 128);

    free(data);
}
