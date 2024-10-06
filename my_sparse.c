#include "spmv.h"


typedef struct {
    int row;   
    int col;   
    double val; 
} SparseMat;

int my_sparse(const unsigned int n, const double mat[],SparseMat sparse[])
{


int cont=0, k=0;


  for (unsigned int i = 0; i < n; i++) {
    sparse[i].row=cont;
      for (unsigned int j = 0; j < n; j++) {
        
          if (mat[i * n + j]!=0){
            cont++;
            sparse[k].col=j;
            sparse[k].val=mat[i * n + j];
            k++;
          }

      }
  }
  sparse[n].row=cont;




  // code your own solver
  return(0);
}
