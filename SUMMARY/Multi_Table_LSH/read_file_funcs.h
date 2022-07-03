#ifndef READ_VECTOR_FILE_FUNCS_H
#define READ_VECTOR_FILE_FUNCS_H

int read_vector_file(char* fPath, int** arr);

int read_vector_file2(char* fPath, int** arr);

int read_hyperplane_tables(int nTables, int nPlanes, int vectorLength, char* fPath,
                           float* hyperplanes);

#endif /* READ_VECTOR_FILE_FUNCS_H */