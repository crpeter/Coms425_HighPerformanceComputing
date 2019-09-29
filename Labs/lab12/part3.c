#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

//void allocate(int** n, int local_n, MPI_Comm comm);

/*void allocate(int** n, int* local_n, int count, int comm_sz, MPI_Comm comm) {
  * n = malloc(*local_n * sizeof(int));
  *local_n = count / comm_sz;
  if (*n == NULL) {
    fprintf(stderr, "Can't allocate vectors\n");
    exit(-1);
  }
  }*/

int main(int argc, char* argv[]) {
  int n, nBuffer, i;
  int my_rank, comm_sz;
  int* local_n_arr, local_n;
  MPI_Comm comm;
  
  MPI_Init(NULL, NULL);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &comm_sz);
  MPI_Comm_rank(comm, &my_rank);
  
  if (my_rank == 0) {
    printf("Enter the number of elts\n");
    scanf("%d", &n);
  }
  
  //allocate(local_n_arr, n / comm_sz, n, comm_sz, comm);

  //Send the value of n just received to other processes
  MPI_Bcast(&n, 1, MPI_INT, 0, comm);
  
  if (my_rank == 0) {
    printf("Process 0\n");
    n = 0;
  } else {
    printf("Process %d: n=%d\n", my_rank, n);
  }


  // Reduce value of n to the min(0, <user entered value>) and send to process 0
  MPI_Reduce(&n, &nBuffer, 1, MPI_INT, MPI_MIN, 0, comm);
  if (my_rank == 0) {
    printf("Process 0 reduced -> n=%d\n\n", n);
  }

  // Reduce value of n to min and send back to all processes
  MPI_Allreduce(&nBuffer, &n, 1, MPI_INT, MPI_MIN, comm);
  
  printf("after all reduce\n");
  int* vec = NULL;
  if (my_rank == 0) {
    vec = malloc(n*sizeof(int));
    for (i = 0; i < n; i++)
      vec[i] = i;
    int a;
    // Scatter value of ver (1 value) to int a.
    MPI_Scatter(vec, local_n, MPI_INT,
		&a, 1, MPI_INT, 0, comm);
    free(vec);
  } else {
    int a;
    MPI_Scatter(vec, 1, MPI_INT,
		&a, 1, MPI_INT, 0, comm);
  }
    
  printf("after scatter\n");
  
  /*double* A = NULL;
  double* help = malloc(sizeof(double));
  help[0] = 1.0;
  if (my_rank == 0) {
    A = malloc(n*sizeof(double));
    MPI_Gather(help, 1, MPI_DOUBLE,
	       A, 1, MPI_DOUBLE, 0, comm);
    for (i = 0; i < n; i++)
      printf("%f ", A[i]);
    printf("\n");
    free(A);
  } else {
    printf("not process 1\n");
    MPI_Send(help, 1, MPI_DOUBLE,
    A, 1, MPI_DOUBLE, 0, comm);
  }
  printf("after gather\n");
  */
  MPI_Finalize();
  return 0;
}
