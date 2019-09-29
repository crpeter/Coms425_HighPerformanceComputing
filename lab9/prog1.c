#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>


/*int main(int argc, char* argv[]) {
  




  return 0;
  }*/


void serSpMatVec(int* I, int* J, double* val, int M, int N, int nz, double* x, double* y) {
  //printf("I: %d, J: %d, val: %d, M: %d", I, J, val, M);
  int i, j, operations = 0;
  
  for (i = 0; i < nz; i++) {
    //printf("y: %f, I: %d, J: %d, val %20.19g\n", y[I[i]], I[i], J[i], val[i]);
    //y[I[i]] = 0.0;
    y[I[i]] += val[i] * x[J[i]];
    //operations++;
    //printf("%f\n", y[I[i]]);
    //for (j = 0; j < N; j++) {
  }
  //printf("\nserSpMat had %d operations\n", operations);
}

void ompSpMatVec(int* I, int* J, double* val, int M, int N, int nz, double* x, double* y, int th) {
  int i, j, operations = 0;
# pragma omp parallel for num_threads(th)
  for (i = 0; i < nz; i++) {
    //printf("y: %f, I: %d, J: %d, val %20.19g\n", y[I[i]], I[i], J[i], val[i]);          
    //y[I[i]] = 0.0;                                                                     
# pragma omp critical
    y[I[i]] += val[i] * x[J[i]];
    //# pragma omp critical    
    //operations++;
  }
  //printf("\nompp had %d operations\n", operations);
}
