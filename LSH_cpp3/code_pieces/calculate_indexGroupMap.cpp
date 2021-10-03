
// Tried doing several tables at a time, but it only made the program slower

void calculate_indexGroupMap(
    // input
    int vectorLength, int numTables,
    int nPoints, float* __restrict__ points,
    int nPlanes, float* __restrict__ hyperplanes,
    int indexGroupMapTableLen,
    // output
    int* __restrict__ indexGroupMap)
{
    const int hyperplanesTableLen = nPlanes * vectorLength;
    #pragma acc data \
      copy(points[numTables*nPoints*vectorLength]) \
      copy(hyperplanes[numTables*nPlanes*vectorLength]) \
      copyout(indexGroupMap[numTables*indexGroupMapTableLen])
    {
        #pragma acc parallel loop collapse(2)
        for (int table = 0; table < numTables; table++) {
            for (int i = 0; i < nPoints; i++) {
                float* hyperplanes2 = hyperplanes + table * hyperplanesTableLen;
                int* indexGroupMap2 = indexGroupMap + table * indexGroupMapTableLen;
                // - Match points -
                // - Calculate hash values, keep track of sizes of each group -
                float* point = &points[i * vectorLength];
                int hashcode = 0;  //  hashcode will be the group index

                // calculate hash value of the i'th point, store resut in indexGroupMap
                #pragma acc loop reduction(|:hashcode)
                for (int j = 0; j < nPlanes; j++) {
                    float* hplane = &hyperplanes2[j * vectorLength];  // first hyperplane
                    // calculate point * hplane
                    float vecMul = 0;
                    #pragma acc loop reduction(+:vecMul)
                    for (int k = 0; k < vectorLength; k++) {
                        vecMul += point[k] * hplane[k];
                    }
                    // set i'th bit to one if point is on "positive" side of hyperplane
                    if (vecMul > 0) {
                        hashcode = hashcode | (1 << j);
                    }
                }
                indexGroupMap2[i] = hashcode;  // save the hashcode
            }
        }
    }
}
