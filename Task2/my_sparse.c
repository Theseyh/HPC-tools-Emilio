#include "spmv.h"

typedef struct {
    int row;   
    int col;   
    double val; 
} SparseMat;

int my_sparse(const unsigned int n, SparseMat sparse[], double vec[], double result[])
{

    
  int j,n1,n2,diff,temp,cont2 = 0;
  for (unsigned int i = 0; i < n; i++){
    n1=sparse[i].row;
    n2=sparse[i+1].row;
    diff=n2-n1;  // the number of element in the row i
    for (j=0;j<diff;j++){
      temp =cont2+j;
      result[i] +=sparse[temp].val * vec[sparse[temp].col];
    }
    cont2+=diff;
  }

    return 0; 
}


