#include "spmv.h"

typedef struct {
    int row;   
    int col;   
    double val; 
} SparseMat;

int my_coo(SparseMat sparse[],const int nnz, double vec[], double result[])
{

  int i,j;  
  double val;
  
// Loop over each non-zero element in the sparse matrix
    for (int k = 0; k < nnz; k++) {
        i = sparse[k].row;      // Row index
        j = sparse[k].col;      // Column index
        val = sparse[k].val; // Non-zero value

        // Perform the multiplication and accumulate the result
        result[i] += val * vec[j];
    }

    return 0; 
}


