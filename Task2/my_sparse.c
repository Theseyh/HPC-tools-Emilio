#include "spmv.h"

typedef struct {
    int row;   // Row index or the starting index of non-zero elements in a row
    int col;   // Column index of the non-zero element
    double val; // Non-zero value
} SparseMat;

int my_sparse(const unsigned int n, SparseMat sparse[], double vec[], double result[])
{
    int j, start, end;

    // Iterate through each row
    for (unsigned int i = 0; i < n; i++) {
        start = sparse[i].row;        // Starting index of non-zero elements in row i
        end = sparse[i + 1].row;      // Ending index (exclusive) for non-zero elements in row i
                                      // (end - start) gives the number of non-zero elements in row i

        // Iterate over the non-zero elements in row i
        for (j = start; j < end; j++) {
            result[i] += sparse[j].val * vec[sparse[j].col]; 
            // Multiply the non-zero value with the corresponding vector element and accumulate in the result
        }
    }

    return 0; 
}
