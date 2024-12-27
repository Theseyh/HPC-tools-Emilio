#include "spmv.h"

int my_dense(const unsigned int n, const double mat[], double vec[], double result[])
{

    // Multiplication matrix-vector
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            result[i] += mat[i * n + j] * vec[j];
        }
    }

    return 0; 
}

