#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


void Allocate_vectors(double** x_pp, double** y_pp, double** z_pp, int n,
                      int local_n, MPI_Comm  comm);
void Get_dims(int* n_p, int* local_n_p,
              int my_rank, int comm_sz, MPI_Comm comm);
void Read_vector(char prompt[], double local_vec[], int n, int local_n,
                 int my_rank, MPI_Comm comm);
void Print_vector(char title[], double local_vec[], int n,
                  int local_n, int my_rank, MPI_Comm comm);
void Vect_mult(double local_A[], double local_x[],
               double local_y[], int n, int local_n,
               MPI_Comm comm);
void Check_for_error(int local_ok, char fname[], char message[],
                     MPI_Comm comm);

int main(int argc, char* argv[]) {
  double* local_A;
  double* local_x;
  double* local_y;
  int n, local_n;
  int my_rank, comm_sz;
  MPI_Comm comm;

  MPI_Init(NULL, NULL);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &comm_sz);
  MPI_Comm_rank(comm, &my_rank);

  Get_dims(&n, &local_n, my_rank, comm_sz, comm);
  Allocate_vectors(&local_A, &local_x, &local_y, n, local_n, comm);

  Read_vector("A", local_A, n, local_n, my_rank, comm);
#  ifdef DEBUG
  Print_vector("A", local_A, n, local_n, my_rank, comm);
#  endif
  
  Read_vector("x", local_x, n, local_n, my_rank, comm);
#  ifdef DEBUG
  Print_vector("x", local_x, n, local_n, my_rank, comm);
#  endif

  Vect_mult(local_A, local_x, local_y, n, local_n, comm);

  Print_vector("y", local_y, n, local_n, my_rank, comm);

  free(local_A);
  free(local_x);
  free(local_y);
  MPI_Finalize();
  return 0;
}


void Allocate_vectors(double**  x_pp  /* out */,
                      double**  y_pp  /* out */,
                      double**  z_pp  /* out */,
                      int       n     /* in  */,
                      int       local_n  /* in  */,
                      MPI_Comm  comm  /* in  */) {
  
  int local_ok = 1;

  *  x_pp = malloc(n*sizeof(double));
  *  y_pp = malloc(n*sizeof(double));
  *  z_pp = malloc(n*sizeof(double));
  
  if (*x_pp == NULL || *y_pp == NULL || *z_pp == NULL) {
    fprintf(stderr, "Can't allocate vectors\n");
    exit(-1);
  }
  Check_for_error(local_ok, "Allocate_arrays",
                  "Can't allocate local arrays", comm);
}


void Get_dims(
              int*      n_p        /* out */,
              int*      local_n_p  /* out */,
              int       my_rank    /* in  */,
              int       comm_sz    /* in  */,
              MPI_Comm  comm       /* in  */) {
  int local_ok = 1;
  
  if (my_rank == 0) {
    printf("Enter the number of elts\n");
    scanf("%d", n_p);
  }
  MPI_Bcast(n_p, 1, MPI_INT, 0, comm);
  if (*n_p <= 0 || *n_p % comm_sz != 0) local_ok = 0;
  Check_for_error(local_ok, "Get_dims",
                  "n must be positive and evenly divisible by comm_sz",
                  comm);
  
  *local_n_p = *n_p/comm_sz;
}


void Read_vector(
                 char      prompt[]     /* in  */,
                 double    local_vec[]  /* out */,
                 int       n            /* in  */,
                 int       local_n      /* in  */,
                 int       my_rank      /* in  */,
                 MPI_Comm  comm         /* in  */) {
  
  double* vec = NULL;
  int i, local_ok = 1;
  
  if (my_rank == 0) {
    vec = malloc(n*sizeof(double));
    if (vec == NULL) local_ok = 0;
    Check_for_error(local_ok, "Read_vector",
                    "Can't allocate temporary vector", comm);
    printf("Enter the vector %s\n", prompt);
    for (i = 0; i < n; i++)
      scanf("%lf", &vec[i]);
    MPI_Scatter(vec, local_n, MPI_DOUBLE,
                local_vec, local_n, MPI_DOUBLE, 0, comm);
    free(vec);
  } else {
    Check_for_error(local_ok, "Read_vector",
                    "Can't allocate temporary vector", comm);
    MPI_Scatter(vec, local_n, MPI_DOUBLE,
                local_vec, local_n, MPI_DOUBLE, 0, comm);
  }
}

void Print_vector(
                  char      title[]     /* in */,
                  double    local_vec[] /* in */,
                  int       n           /* in */,
                  int       local_n     /* in */,
                  int       my_rank     /* in */,
                  MPI_Comm  comm        /* in */) {
  double* vec = NULL;
  int i, local_ok = 1;
  
  if (my_rank == 0) {
    vec = malloc(n*sizeof(double));
    if (vec == NULL) local_ok = 0;
    Check_for_error(local_ok, "Print_vector",
                    "Can't allocate temporary vector", comm);
    MPI_Gather(local_vec, local_n, MPI_DOUBLE,
               vec, local_n, MPI_DOUBLE, 0, comm);
    printf("\nThe vector %s\n", title);
    for (i = 0; i < n; i++)
      printf("%f ", vec[i]);
    printf("\n");
    free(vec);
  }  else {
    Check_for_error(local_ok, "Print_vector",
                    "Can't allocate temporary vector", comm);
    MPI_Gather(local_vec, local_n, MPI_DOUBLE,
               vec, local_n, MPI_DOUBLE, 0, comm);
  }
}

void Vect_mult(
               double    local_A[]  /* in  */,
               double    local_x[]  /* in  */,
               double    local_y[]  /* out */,
               int       n          /* in  */,
               int       local_n    /* in  */,
               MPI_Comm  comm       /* in  */) {
  double* x;
  int local_i;
  int local_ok = 1;
  
  x = malloc(n*sizeof(double));
  if (x == NULL) local_ok = 0;
  Check_for_error(local_ok, "Mat_vect_mult",
                  "Can't allocate temporary vector", comm);
  MPI_Allgather(local_x, local_n, MPI_DOUBLE,
                x, local_n, MPI_DOUBLE, comm);
  
  for (local_i = 0; local_i < local_n; local_i++) {
    local_y[local_i] = local_A[local_i] * local_x[local_i];
  }
  free(x);
}


// Helper



/*-------------------------------------------------------------------
 * Function:  Check_for_error
 * Purpose:   Check whether any process has found an error.  If so,
 *            print message and terminate all processes.  Otherwise,
 *            continue execution.
 * In args:   local_ok:  1 if calling process has found an error, 0
 *               otherwise
 *            fname:     name of function calling Check_for_error
 *            message:   message to print if there's an error
 *            comm:      communicator containing processes calling
 *                       Check_for_error:  should be MPI_COMM_WORLD.
 *
 * Note:
 *    The communicator containing the processes calling Check_for_error
 *    should be MPI_COMM_WORLD.
 */
void Check_for_error(
                     int       local_ok   /* in */,
                     char      fname[]    /* in */,
                     char      message[]  /* in */,
                     MPI_Comm  comm       /* in */) {
  int ok;
  
  MPI_Allreduce(&local_ok, &ok, 1, MPI_INT, MPI_MIN, comm);
  if (ok == 0) {
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);
    if (my_rank == 0) {
      fprintf(stderr, "Proc %d > In %s, %s\n", my_rank, fname,
              message);
      fflush(stderr);
    }
    MPI_Finalize();
    exit(-1);
  }
}
