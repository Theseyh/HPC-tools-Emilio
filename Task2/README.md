In this folder, we have the script "Lebatsh.sh", which is used for running batches of sparse matrix-vector multiplication (.spmv) operations compiled with different optimization levels. The folder contains separate directories for batches compiled with ICC and GCC.

The **.spmv** format executes all the functions: **my_coo**, **my_csr**, **my_csc**, and the reference functions from **GSL** and **MKL**.


The two tables present the performance of sparse matrix operations (my_dense, my_coo, my_csr, and my_csc) for 16384 × 16384 matrices with 10% non-zero elements, evaluated using different optimization levels for GCC (Table 1) and ICC (Table 2). All values are given in milliseconds (ms), and the data shown was derived from an average of five measurements for each combination of matrix format and optimization level.

 Table 1: GCC Benchmark - 16384 × 16384 matrices, 10% non-zero elements 
|   in ms    | O0     | O2-novec | O3-vec | Ofast-vec | Ref |
|-------------|--------|----------|--------|-----------|-----|
| my_dense    |    401.4    |    401.0      |   401.0     |      401.0     |   129.0  |
| my_coo      |    88.8    |      88.0      |    89.0    |    89.0     |  87.0   |
| my_csr      |    42.6    |       42.8     |     42.8   |   42.8    |   38.6  |
| my_csc      |     40.4   |      40.4      |    40.2    |    40.4     |   36.0  |




Table 2: ICC Benchmark - 16384 × 16384 matrices, 10% non-zero elements 
|     in ms        | O0     | O2-novec | O3-vec | fast | Ref |
|-------------|--------|----------|--------|-----------|-----|
| my_dense    |    902.8    |    399.8      |   149.2    |      158.6     |   129.0  |
| my_coo      |    147.8    |      89.4      |    89.4    |    89.2     |  87.0   |
| my_csr      |    98.0    |       42.4     |     42.4   |   35.0   |   38.6  |
| my_csc      |     105.6   |      41.0      |    41.0    |    39.0    |   36.0  |

### Comments on the Benchmark Tables



**Table 1: GCC Benchmark**
- The GCC benchmark highlights the performance of four different matrix formats across various optimization levels (O0, O2-novec, O3-vec, Ofast-vec).
- **My_dense** consistently performed similarly across the different optimization levels, with an average value of around 401 ms, significantly higher than the reference (129 ms).
- **My_coo** showed a slight variation with optimization, ranging from 88.0 ms (O2-novec) to 89.0 ms (O3-vec and Ofast-vec), but remained quite close to the reference (87.0 ms).
- **My_csr** demonstrated consistent results across optimizations, with values around 42.6 ms to 42.8 ms, close to the reference of 38.6 ms.
- **My_csc** had a stable performance with optimization levels O0 and Ofast-vec at 40.4 ms, and a slight decrease to 40.2 ms at O3-vec. The reference value was 36.0 ms.

**Table 2: ICC Benchmark**
- The ICC benchmark shows a more varied set of results, with the **my_dense** computation being notably slower than in GCC, ranging from 902.8 ms at O0 to 149.2 ms at O3-vec and 158.6 ms at fast, compared to the reference value of 129 ms.
- The **my_coo** format performed relatively stable, with values ranging from 147.8 ms (O0) to 89.4 ms (O3-vec, fast), slightly improving with optimization.
- **My_csr** computation also showed improvement with optimizations, from 98.0 ms at O0 to 42.4 ms at O2-novec and O3-vec, and down to 35.0 ms at the "fast" setting, which is below the reference of 38.6 ms.
- **My_csc** also showed a good improving performance, with differences between optimization levels, ranging from 105.6 ms (O0) to 39.0 ms at the "fast" setting. The reference was again 36.0 ms.

