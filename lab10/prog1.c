#include <stdio.h>
#include "omp.h"

void serMM2CSR(int* I, int* J, double* valMM, int M, int N, int nz, double* valCSR, int* col_ind, int* row_ptr);

void serSpMatVec(double* val, int* col_ind, int* row_ptr, int M, int N, int nz, double* x, double* y);

void ompSpMatVec(double* val, int* col_ind, int* row_ptr, int M, int N, int nz, double* x, double* y, int th);

void serMM2CSR(int* I, int* J, double* valMM, int M, int N, int nz, double* valCSR, int* col_ind, int* row_ptr) {
  int i, count = 0, row_count = 1, last_col = 0;
  
  row_ptr[0] = 0;
  for (i = 0; i < nz; i++) {
    //printf("I: %d, J: %d, val %20.19g\n", I[i], J[i], valMM[i]);
    valCSR[i] = valMM[i];
    col_ind[i] = J[i];
    if (I[i] > I[i-1]) {
      row_ptr[row_count] = i;
      row_count++;
    } else if (i == nz-1) {
      row_ptr[row_count] = i+1;
    }
  }
//  printf("###\n\n\n###\n");
//  for (i = 0; i < M+1; i++) {
//    printf("%d: %d\n", i, row_ptr[i]);
//  }
  /*printf("\n");
  for (i = 0; i < 4; i++) {
    printf("%d: %f\n", i, valCSR[i]);
    }*/
}


void serSpMatVec(double* val, int* col_ind, int* row_ptr, int M, int N, int nz, double* x, double* y) {
  int i, j, row_count = 0;
  
  for (i = 0; i < N; i++) y[i] = 0.0;
  
  for (i = 0; i < N; i++) {
    for (j = row_ptr[i]; j < row_ptr[i+1]; j++) {
      y[i] += val[j] * x[j];
    }
  }
  /*
  printf("###\n\n\n###\n");
  for (i = 0; i < 10; i++) {
    printf("%d: %f\n", i, y[i]);
    }*/
}


void ompSpMatVec(double* val, int* col_ind, int* row_ptr, int M, int N, int nz, double* x, double* y, int th) {
  int i, j, row_count = 0;
  
  for (i = 0; i < N; i++) y[i] = 0.0;
  
# pragma omp parallel num_threads(th) private(i) shared(N, row_ptr, y, val, x)
  for (i = 0; i < N; i++) {
    //# pragma omp for
    for (j = row_ptr[i]; j < row_ptr[i+1]; j++) {
      y[i] += val[j] * x[j];
    }
  }
  /*
  printf("###\n\n\n###\n");
  for (i = 0; i < 10; i++) {
    printf("%d: %f\n", i, y[i]);
  }
  */
}
