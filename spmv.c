#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>      // CBLAS in GSL (the GNU Scientific Library)
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
#include "timer.h"
#include "spmv.h"

#define DEFAULT_SIZE 1024
#define DEFAULT_DENSITY 0.25


typedef struct {
    int row;   
    int col;   
    double val; 
} SparseMat;


unsigned int populate_sparse_matrix(double mat[], unsigned int n, double density, unsigned int seed)
{
  unsigned int nnz = 0;

  srand(seed);

  for (unsigned int i = 0; i < n * n; i++) {
    if ((rand() % 100) / 100.0 < density) {
      // Get a pseudorandom value between -9.99 e 9.99
      mat[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
      nnz++;
    } else {
      mat[i] = 0;
    }
  }

  return nnz;
}

unsigned int populate_vector(double vec[], unsigned int size, unsigned int seed)
{
  srand(seed);

  for (unsigned int i = 0; i < size; i++) {
    vec[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
  }

  return size;
}

int is_nearly_equal(double x, double y)
{
  const double epsilon = 1e-5 /* some small number */;
  return fabs(x - y) <= epsilon * fabs(x);
  // see Knuth section 4.2.2 pages 217-218
}

unsigned int check_result(double ref[], double result[], unsigned int size)
{
  for(unsigned int i = 0; i < size; i++) {
    if (!is_nearly_equal(ref[i], result[i]))
      return 0;
  }

  return 1;
}

int main(int argc, char *argv[])
{
  int size;        // number of rows and cols (size x size matrix)
  double density;  // aprox. ratio of non-zero values

  if (argc < 2) {
    size = DEFAULT_SIZE;
    density = DEFAULT_DENSITY;
  } else if (argc < 3) {
    size = atoi(argv[1]);
    density = DEFAULT_DENSITY;
  } else {
    size = atoi(argv[1]);
    density = atoi(argv[2]);
  }

  double *mat, *vec, *refsol, *mysol;

  mat = (double *) malloc(size * size * sizeof(double));
  vec = (double *) malloc(size * sizeof(double));
  refsol = (double *) malloc(size * sizeof(double));
  mysol = (double *) malloc(size * sizeof(double));

  unsigned int nnz = populate_sparse_matrix(mat, size, density, 1);
  populate_vector(vec, size, 2);

  printf("Matriz size: %d x %d (%d elements)\n", size, size, size*size);
  printf("%d non-zero elements (%.2lf%%)\n\n", nnz, (double) nnz / (size*size) * 100.0);

  //
  // Dense computation using CBLAS (eg. GSL's CBLAS implementation)
  //
  printf("Dense computation\n----------------\n");

  timeinfo start, now;
  timestamp(&start);

  cblas_dgemv(CblasRowMajor, CblasNoTrans, size, size, 1.0, mat, size, vec, 1, 0.0, refsol, 1);

  timestamp(&now);
  printf("Time taken by CBLAS dense computation: %ld ms\n", diff_milli(&start, &now));

  //
  // Using your own dense implementation
  //
  timestamp(&start);

  my_dense(size, mat, vec, mysol);

  timestamp(&now);
  printf("Time taken by my dense matrix-vector product: %ld ms\n", diff_milli(&start, &now));

  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");


  //
  // Let's try now SpMV: Sparse Matrix - Dense Vector computation
  //

  // Convert mat to a sparse format: CSR
  // Use the gsl_spmatrix struct as datatype

  //
  // Sparse computation using GSL's sparse algebra functions
  //

    //
  // Convert mat to a sparse format: CSR
  // Use the gsl_spmatrix struct as datatype
  //

  gsl_spmatrix *sparse_mat = gsl_spmatrix_alloc(size, size);
  
  // Populate the sparse matrix from the dense matrix
  for ( int i = 0; i < size; i++) {
    for ( int j = 0; j < size; j++) {
      double value = mat[i * size + j];
      if (value != 0.0) {
        gsl_spmatrix_set(sparse_mat, i, j, value);
      }
    }
  }



  //
  // Sparse computation using GSL's sparse algebra functions
  //
  printf("\nSparse computation\n----------------\n");

  gsl_vector *gsl_vec = gsl_vector_alloc(size);
  gsl_vector *gsl_result = gsl_vector_alloc(size);

  // Copy data from the dense vector to the GSL vector
  for ( int i = 0; i < size; i++) {
    gsl_vector_set(gsl_vec, i, vec[i]);
  }

  // Sparse matrix-vector multiplication using GSL's sparse algebra functions
  timestamp(&start);

  gsl_spblas_dgemv(CblasNoTrans, 1.0, sparse_mat, gsl_vec, 0.0, gsl_result);

  timestamp(&now);
  printf("Time taken by GSL sparse matrix-vector product: %ld ms\n", diff_milli(&start, &now));

  // Copy the result back to the mysol array for comparison
  for ( int i = 0; i < size; i++) {
    mysol[i] = gsl_vector_get(gsl_result, i);
    //printf ("%lf et ref %lf\n",mysol[i], refsol[i]);
  }
  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");
  printf("test");
  // Free GSL resources
  gsl_spmatrix_free(sparse_mat);
  gsl_vector_free(gsl_vec);
  gsl_vector_free(gsl_result);


  //
  // Your own sparse implementation
  //
  printf("test2");
  SparseMat *sparse = (SparseMat *)malloc(size * size * sizeof(SparseMat));

  // Appel de la fonction
  int non_zero_elements = my_sparse(size, mat, sparse);
  printf("%d",non_zero_elements);

  // Affichage des éléments non nuls
  for (unsigned int i = 0; i < nnz; i++) {// pas nnz mais sparse[i].row
    //printf("Element non nul: row=%d, col=%d, val=%lf\n", sparse[i].row, sparse[i].col, sparse[i].val);
  }
  int k1,n1,nplusun,diff,cont2 = 0;
   int j;
  printf("\n%d\t%d\n",sparse[size-1].row,sparse[nnz-2].col);
  for ( int i = 0; i < size; i++){
    n1=sparse[i].row;
    nplusun=sparse[i+1].row;
    //printf("\n%d\n",diff);
    diff=nplusun-n1;
    //printf("\n%d\n",diff);
    for (j=0;j<diff;j++){
    printf("Element non nul: row=%d, col=%d, val=%lf\n", i, sparse[cont2+j].col, sparse[cont2+j].val);
    //printf("nnz = %d\n",nnz);
    printf("cont2 = %d\n",cont2);
    k1++;
    }
    cont2+=diff;
    printf("\n%d\t%d\n",k1,nnz);
  }
  

  // Compare times (and computation correctness!)


  // Free resources
  free(mat);
  free(vec);
  free(refsol);
  free(mysol);

  return 0;
}
