

# SpMV: Sparse Matrix-Vector Product

This document describes the code implementation and benchmarks for performing Sparse Matrix-Vector product (SpMV) in both dense and sparse formats.

## Libraries Used
The code uses the following libraries:
- **GSL (GNU Scientific Library)** for sparse and dense matrix-vector operations.
  - Dense matrix-vector product: `cblas_dgemv()`, requires linking against `libgslcblas`.
  - Sparse matrix-vector product: `gsl_spblas_dgemv()`, requires linking against `libgsl`.

Other CBLAS implementations can be used for the dense product, such as `-lopenblas` instead of `-lgslcblas`.

### Execution Instructions

To execute the code, follow these steps:

1. **Build the code**:
   ```bash
   make
   ```

2. **Run the program**:
   ```bash
   ./spmv
   ```

### Overview

In this work, we implemented four different methods to perform matrix-vector multiplication:

1. **CBLAS Dense**: 
   Utilizes the `cblas_dgemv()` function from the CBLAS library to multiply a dense matrix by a vector.

2. **My Dense**: 
   A custom implementation of dense matrix-vector multiplication.

3. **GSL Sparse**: 
   Uses the GSL library’s sparse matrix-vector multiplication function, `gsl_spblas_dgemv()`, after converting a dense matrix into sparse format.

4. **My Sparse**: 
   A custom implementation of sparse matrix-vector multiplication using the CSR (Compressed Sparse Row) format.



## Performance Comparison

The table below summarizes the computation times for various methods, based on 10 executions of the program. The matrix used was of size 1024x1024, with a vector of size 1024. Additionally, 1/4 of the matrix's values were non-null, making it a sparse matrix for the sparse methods.

| Computation Type        | Mean Time (μs) | Variance (μs²) | Median Time (μs) |
|-------------------------|----------------|----------------|------------------|
| **CBLAS Dense**          | 1223.2         | 402,562.76     | 793.5            |
| **My Dense**             | 1842.1         | 258,252.49     | 1617.0           |
| **GSL Sparse**           | 299.5          | 2,965.05       | 289.5            |
| **My Sparse**            | 403.5          | 5,743.05       | 371.0            |

### Analysis

- **My Dense** computations are the slowest, which is expected as the code is not optimized for dense matrix-vector multiplication. This leads to inefficient resource use and slower execution times.

- **CBLAS Dense** computations are faster than my dense implementation. CBLAS provides an optimized method for performing dense matrix-vector multiplication, which leads to better performance, but it still has to deal with the overhead of processing a dense matrix.

- **My Sparse** computations come next in terms of speed. This method uses a sparse matrix in CSR (Compressed Sparse Row) format for matrix-vector multiplication. While it benefits from the sparse format, the implementation isn't fully optimized for memory access, which impacts performance compared to the GSL method.

- **GSL Sparse** computations are the fastest. GSL’s optimized sparse matrix-vector multiplication leverages efficient memory access patterns and algorithms to fully exploit the sparse matrix structure, leading to the best performance.

This logical progression shows that more optimized methods and the use of sparse matrix formats lead to faster computations, with GSL being the most efficient for sparse matrix operations.





