#include "spmv.h"
#include <stdio.h>

typedef struct {
    int row;   // Row index
    int col;   // Starting index of each column
    double val; // Non-zero value
} SparseMat;

int my_csc(const unsigned int n, SparseMat sparse[], double vec[], double result[])
{

    // Iterate through each column
    for (unsigned int j = 0; j < n; j++) {
        // Retrieve the range of indices for this column
        int start = sparse[j].col;
        int end = sparse[j+1].col;

        // Iterate over the non-zero elements of this column
        for (int k = start; k < end; k++) {
            // Accumulate in the result for the corresponding row
            result[sparse[k].row] += sparse[k].val * vec[j];

        }
    }

    return 0; 
}
