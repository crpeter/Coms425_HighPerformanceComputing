/* 
*   Matrix Market I/O example program
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
*   and copies it to stdout.  This porgram does nothing useful, but
*   illustrates common usage of the Matrix Matrix I/O routines.
*   (See http://math.nist.gov/MatrixMarket for details.)
*
*   Usage:  a.out [filename] > output
*/

#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "prog1.h"
#include "omp.h"
int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   
    int i, *I, *J;
    double *val;

    if (argc < 2) {
      fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
    } else { 
        if ((f = fopen(argv[1], "r")) == NULL) 
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
    //int num_threads = strtol(argv[2], NULL, 10);


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reseve memory for matrices */

    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    //void serSpMatVec(int* I, int* J, double* val, int M, int N, int nz, double* x, double* y)
    double* x = malloc(M * sizeof(double));
    double* y = malloc(M * sizeof(double));
    for(i = 0; i < M; i++) {
      x[i] = (i % 100) * ((i-1) * (i+1) - (i-2) / (i+2));
      y[i] = 0.0;
    }
    printf("M: %d, N: %d, nz: %d\n", M, N, nz);
    double totalTime;
    //for(i=0; i < 16; i++) {
    totalTime = omp_get_wtime();
    //printf("time: %d\n", omp_get_wtime());
    serSpMatVec(I, J, val, M, N, nz, x, y);
    //printf("now time: %d\n", omp_get_wtime());
    printf("Cereal time: %f\n", omp_get_wtime() -totalTime);
    //}
    for(i = 0; i < M; i++) {
      ;//printf("%5g ", y[i]);
    }
    printf("\n");
    for(i = 0; i < M; i++) {
      y[i] = 0.0;
    }
    for(i=1; i < 17; i++) {
    totalTime = omp_get_wtime();
    ompSpMatVec(I, J, val, M, N, nz, x, y, i);
    printf("%3d Threads, Time: %f\n", i, omp_get_wtime() -totalTime);
    }
    //printf("\n###\n###omp ones###\n###\n");
    for(i = 0; i < M; i++) {
      ;//printf("%5g ", y[i]);
    }
    printf("\n");
    /*mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    
    for (i=0; i<nz; i++)
        fprintf(stdout, "%d %d %20.19g\n", I[i]+1, J[i]+1, val[i]);
    */
    return 0;
}

