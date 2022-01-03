#include <iostream>


int main(int argc, char** argv)
{
    if (argc < 3) {
        fprintf(stderr, "Program needs to be given two files as arguments:");
        fprintf(stderr, " Binary file with feature descriptors");
        fprintf(stderr, " and a binary file with the ground truth table\n");
        return -1;
    }

    float* points;
    int* dimensions;
    // int* buffer;

    char* fPath = argv[1];
    printf("%s\n", fPath);

    // Opening file containing feature vectors

    FILE* fp = fopen(fPath, "rb");

    if (fp == NULL) {
        fprintf(stderr, "Error opening %s for reading\n", fPath);
        return -1;
    }

    // Get the size of the feature vectors
    fseek(fp, 0L, SEEK_END);
    int filesize = ftell(fp);
    rewind(fp);

    // buffer = (int*) malloc(filesize);

    // Allocates space for the point
    points = (float*) malloc(filesize - (filesize / 129));
    dimensions = (int*) malloc(filesize / 129);

    int nPoints = (filesize / 129) / 4;

    int ret = 0;
    for (int i = 0; i < nPoints; i++) {
        ret = fread(dimensions+i, 4, 1, fp);
        if (ret != 1) {
            fprintf(stderr, "Premature end of file\n");
            fprintf(stderr, " Error when reading point nr: %d.\n", i);
            fprintf(stderr, " Only read %d elements should have read 1\n", ret);
            break;
        }
        ret = fread(points+128*i, 4, 128, fp);
        if (ret != 128) {
            fprintf(stderr, "Premature end of file\n");
            fprintf(stderr, " Error when reading point nr: %d.\n", i);
            fprintf(stderr, " Only read %d elements should have read 128\n", ret);
            break;
        }
    }

    fclose(fp);


    // Opening file containing ground truth table


    fp = fopen(argv[2], "rb");

    // Get the size of the feature vectors
    fseek(fp, 0L, SEEK_END);
    filesize = ftell(fp);
    rewind(fp);
    
    int nEntries = (filesize - filesize / 101) / 4;

    int* dimensions_gt = (int*) malloc((nEntries / 100) * 4);
    int* ground_truth = (int*) malloc(nEntries * 4);

    // printf("nEntries: %d\n", nEntries);
    printf("Dim entries: %d\n", (nEntries / 100));
    printf("Match entries: %d\n", nEntries);

    // ret = fread(dimensions_gt, 4, 1, fp);
    // ret = fread(dimensions_gt, 4, 100, fp);
    // ret = fread(dimensions_gt, 4, 1, fp);
    // printf("dimensions = %d\n", dimensions_gt[0]);

    printf("File size (bytes): %d\n", filesize);
    printf("File size (words): %d\n", filesize/4);

    // for (int i = 0; i < 101; i++) {
    //     int tmp = 0;
    //     fread(&tmp, 4, 1, fp);
    //     printf("%d\n", tmp);
    // }

    ret = 0;
    for (int i = 0; i < nEntries / 101; i++) {
        ret = fread(dimensions_gt + i, 4, 1, fp);
        if (ret != 1) {
            fprintf(stderr, "Premature end of file '%s'\n", argv[2]);
            fprintf(stderr, " Error when reading part nr: %d.\n", i);
            fprintf(stderr, " Only read %d elements should have read 1\n", ret);
            break;
        }
        ret = fread(ground_truth + dimensions_gt[i]*i, 4, dimensions_gt[i], fp);
        if (ret != dimensions_gt[i]) {
            fprintf(stderr, "Premature end of file '%s'\n", argv[2]);
            fprintf(stderr, " Error when reading part nr: %d.\n", i);
            fprintf(stderr, " Only read %d elements should have read 100\n", ret);
            break;
        }
    }

    fclose(fp);


    // // // --- Print some examples values found in the feature descriptor input ---

    // for (int i = 0; i < 20; i++) {
    //     printf("%d ", dimensions[i]);
    // }
    // printf("\n");

    printf("\nSome points");
    printf("------------------\n");
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 128; j++) {
            printf("%f\n", points[i*128 + j]);
        }
        printf("------------------\n");
    }

    // // --- Print maximum value found for the feature descriptor inputs ---
    // int my_max = 0;
    // for (int i = 0; i < nPoints*128; i++) {
    //     my_max = (points[i] > my_max) ? points[i] : my_max;
    // }
    // printf("%d\n", my_max);

    free(points);
    free(dimensions);
    free(dimensions_gt);
    free(ground_truth);
}
