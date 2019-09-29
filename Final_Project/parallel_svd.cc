#include "mpi.h"
#include <cblas.h>
#include <iostream>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include "omp.h"
#define TALLOC(n,typ) (typ*)malloc(n*sizeof(typ))

//#define MATRIX_SIZE 5
#define MAX_PROCESSES 200
#define EL_Type double
#define MATRIX_SIZE 100

using namespace std;

//Function to compute qr - decomposition of matrix
void computeQR(bool verbose, bool definition, float **matrixR, float **matrixQ, int rank, int lines);

//Parallel matrix multiply by vector
void Mat_vect_mult(double local_A[], double local_x[],double local_y[], int local_m, int n, int local_n,MPI_Comm comm);

typedef struct sMtx {
    int     dim_x, dim_y;
    EL_Type *m_stor;
    EL_Type **mtx;
} *Matrix, sMatrix;
 
typedef struct sRvec {
    int     dim_x;
    EL_Type *m_stor;
} *RowVec, sRowVec;

//Create a new Matrix 
Matrix NewMatrix( int x_dim, int y_dim )
{
    int n;
    Matrix m;
    m = (sMtx*)TALLOC( 1, sMatrix);
    n = x_dim * y_dim;
    m->dim_x = x_dim;
    m->dim_y = y_dim;
    m->m_stor = (EL_Type*)TALLOC(n, EL_Type);
    m->mtx = (EL_Type**)TALLOC(m->dim_y, EL_Type *);
    for(n=0; n<y_dim; n++) {
        m->mtx[n] = m->m_stor+n*x_dim;
    }
    return m;
}

//Set row of a matrix
void MtxSetRow(Matrix m, int irow, EL_Type *v)
{
    int ix;
    EL_Type *mr;
    mr = m->mtx[irow];
    for(ix=0; ix<m->dim_x; ix++)
        mr[ix] = v[ix];
}

//Init a Matrix
Matrix InitMatrix( int x_dim, int y_dim, EL_Type **v)
{
    Matrix m;
    int iy;
    m = NewMatrix(x_dim, y_dim);
    for (iy=0; iy<y_dim; iy++) 
        MtxSetRow(m, iy, v[iy]);
    return m;
}

//Display a Matrix
void MtxDisplay( Matrix m )
{
    int iy, ix;
    const char *sc;
    for (iy=0; iy<m->dim_y; iy++) {
        printf("   ");
        sc = " ";
        for (ix=0; ix<m->dim_x; ix++) {
            printf("%s %f", sc, m->mtx[iy][ix]);
            sc = ",";
        }
        printf("\n");
    }
    printf("\n");
}

//Matrix multiply and add rows
void MtxMulAndAddRows(Matrix m, int ixrdest, int ixrsrc, EL_Type mplr)
{
    int ix;
    EL_Type *drow, *srow;
    drow = m->mtx[ixrdest];
    srow = m->mtx[ixrsrc];
    for (ix=0; ix<m->dim_x; ix++) 
        drow[ix] += mplr * srow[ix];
//	printf("Mul row %d by %d and add to row %d\n", ixrsrc, mplr, ixrdest);
//	//	MtxDisplay(m);
}

//Matrix mult by vec
void MtxMulVec(Matrix m, EL_Type *vec, EL_Type *out) {
  int i, j, dim_y = m->dim_y-1;
#pragma omp parallel for
  for (i = 0; i < dim_y; i++) {
    out[i] = 0.0;
    for (j = 0; j < m->dim_x; j++) {
      out[i] += m->mtx[i][j] + vec[i];
    }
  }
}

void MtxSwapRows( Matrix m, int rix1, int rix2)
{
    EL_Type *r1, *r2, temp;
    int ix;
    if (rix1 == rix2) return;
    r1 = m->mtx[rix1];
    r2 = m->mtx[rix2];
    for (ix=0; ix<m->dim_x; ix++)
        temp = r1[ix]; r1[ix]=r2[ix]; r2[ix]=temp;
//	printf("Swap rows %d and %d\n", rix1, rix2);
//	//	MtxDisplay(m);
}

//Normalize row of matrix
void MtxNormalizeRow( Matrix m, int rix, int lead)
{
    int ix;
    EL_Type *drow;
    EL_Type lv;
    drow = m->mtx[rix];
    lv = drow[lead];
    for (ix=0; ix<m->dim_x; ix++)
        drow[ix] /= lv;
    //printf("Normalize row %d\n", rix);
		//MtxDisplay(m);
}


#define MtxGet( m, rix, cix ) m->mtx[rix][cix]

//Reduce Matrix to RE form
void MtxToReducedREForm(Matrix m)
{
    int lead;
    int rix, iix;
    EL_Type lv;
    int rowCount = m->dim_y;
 
    lead = 0;
    for (rix=0; rix<rowCount; rix++) {
        if (lead >= m->dim_x)
            return;
        iix = rix;
        while (0 == MtxGet(m, iix,lead)) {
            iix++;
            if (iix == rowCount) {
                iix = rix;
                lead++;
                if (lead == m->dim_x)
                    return;
            }
        }
        MtxSwapRows(m, iix, rix );
        MtxNormalizeRow(m, rix, lead );
        for (iix=0; iix<rowCount; iix++) {
            if ( iix != rix ) {
                lv = MtxGet(m, iix, lead );
                MtxMulAndAddRows(m,iix, rix, -lv) ;
            }
        }
        lead++;
    }
}

//Compute eigenvectors of matrix
//Results returned in vecBuff param
void EigenVectorNormalize(Matrix m, EL_Type* vecBuff, EL_Type eigenValue){
  int i, j, dim_y = m->dim_y;
	EL_Type sum_of_squares = 0;
#pragma omp parallel for
	for(i = 0; i < dim_y; i++){
		//printf("Augmented Value: %f\n",m->mtx[i][m->dim_x - 1]);
		
		EL_Type check = 1.0/m->mtx[i][m->dim_x - 1];
		if(check < 0){
			check = check * -1.0;
		}
		
		if(check == eigenValue || m->mtx[i][m->dim_x - 1] == 1.0){
			
			
			continue;
		}
		sum_of_squares += m->mtx[i][m->dim_x - 1] * m->mtx[i][m->dim_x - 1];
	}
	EL_Type normalization_scalar = 0.0;
	if(sum_of_squares == 0.0){
		normalization_scalar = 0.0;
	}
	else{
		normalization_scalar = 1.0/(sqrt(sum_of_squares));
	}
	//printf("Norm Scalar: %f\n", normalization_scalar);
#pragma omp parallel for
	for(j = 0; j < dim_y; j++){
		EL_Type check2 = 1.0/(m->mtx[j][m->dim_x - 1]);
		//printf("check2: %f    %f\n", check2 * check2, eigenValue * eigenValue);
		
		if((((check2 * check2) - (eigenValue * eigenValue)) < .000001 && check2*check2 - eigenValue*eigenValue > -.000001) || m->mtx[j][m->dim_x-1] == 1.0 || m->mtx[j][m->dim_x-1] == -1){
		//	printf("hereeeee\n");
			vecBuff[j] = 0.0;
		}
		else{
			vecBuff[j] = normalization_scalar * m->mtx[j][m->dim_x - 1];
		}
	
	}

}

//Print a vector
void printVec(EL_Type* vec, int n){
	int i;
	for(i = 0; i < n; i++){
		printf("Index %d: %f\n", i, vec[i]);
	}

}

//Subtract an Eigen value from the diaganols of a Matrix
void subtractEigenValue(Matrix m, EL_Type eigenValue){
	int i, j;
	for(i = 0; i < m->dim_y; i++){
		EL_Type cur = m->mtx[i][i];
		m->mtx[i][i] = cur - eigenValue;
	}
}

//Create a matrix as a float* pointer (for getting eigenvalues in parallel)
float** create_matrix( int numrows, int numcols){
    float *buffer=new float[numrows*numcols];
    float **data=new float*[numrows];
    for(int i=0;i<numrows;++i) data[i]=buffer+i*numcols;
    
    return *&data;
}


int main(int argc, char* argv[])
{

    int my_rank, comm_sz, i,j,k,l,m,lines;
    float **matrixA,**matrixR,**matrixQ, **matrixOriginal, *eigenValues;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0) printf("Matrix size: %d\n", MATRIX_SIZE);
    srand((unsigned)time(NULL));

    Matrix m1;

    static EL_Type** im = TALLOC(MATRIX_SIZE, EL_Type*);

    int p, q;
    for(p = 0; p < MATRIX_SIZE; p++){
	EL_Type* row = TALLOC(MATRIX_SIZE, EL_Type);
	for(q = 0; q < MATRIX_SIZE; q++){
		int rand_num = rand() % 100;
		row[q] = (EL_Type)rand_num;
	}
	im[p] = row;
    }

   
    m1 = InitMatrix( MATRIX_SIZE + 1, MATRIX_SIZE, im );

    int z;
    for(z = 0; z < m1->dim_y; z++){
	m1->mtx[z][MATRIX_SIZE] = 1.0;

    }    

    EL_Type recvBuff;
    
    //Bcast lines in
    lines = MATRIX_SIZE;
    MPI_Bcast((void *)&lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
    EL_Type* all_eigenvalues = (EL_Type*)TALLOC(lines, EL_Type);

    matrixQ = create_matrix(lines,lines);
    matrixR = create_matrix(lines,lines);
    matrixA = create_matrix(lines,lines);
    //if(my_rank == 0) {
        //Create Matrix R
        for(i = 0; i < lines; i++){
            for(j = 0; j < lines; j++){
                matrixR[i][j] = m1->mtx[i][j];
            }
        }
        for(i = 0; i < lines; i++){
            for(j = 0; j < lines; j++){
                matrixA[i][j] = matrixR[i][j];
            }
        }
    //}
    
    
    int a = 20;
    int b = (int)ceil(sqrt( lines ));
    double totalTime;
    totalTime = omp_get_wtime();
    for(i = 0; i < std::min(a,b); i++) {
      for(j = 0; j < lines; j++) {
        for(k = 0; k < lines; k++) {
          matrixR[j][k] = matrixA[j][k];
          matrixQ[j][k] = 0.0;
        }
      }
      computeQR(false, false, matrixR, matrixQ, my_rank, lines);
      // A = Q'*R
      for(k = 0; k < lines; k++){
        for(l = 0; l < lines; l++){
          float tm = 0;
          for(m = 0; m < lines; m++){
            tm += matrixQ[m][k]*matrixA[m][l];
          }
          matrixR[k][l] = tm;
        }
      }
      // R = A*Q
      for(k = 0; k < lines; k++){
        for(l = 0; l < lines; l++){
          float tm = 0;
          for(m = 0; m < lines; m++){
            tm += matrixR[k][m]*matrixQ[m][l];
          }
          matrixA[k][l] = tm;
        }
      }
    }
    if(my_rank == 0) {
      bool zero = false;
      int values = 0;
      eigenValues = (float*)malloc(lines);
      for(k = 0; k < lines; k++) {
          /**
           *
           * humanCentipad wont read 0.0
           *
           *
           */
        if (matrixA[k][k] == 0.0 && !zero) {
          zero = true;
          all_eigenvalues[k] = matrixA[k][k];
          values++;
        } else {
          all_eigenvalues[k] = matrixA[k][k];
          values++;
        }
      }
      
      for (k = 0; k < values; k++) {
        //printf("eig %d: %20f\n", k, all_eigenvalues[k]);
      }
    }
   
    MPI_Scatter(all_eigenvalues, 1, MPI_DOUBLE, &recvBuff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //printf("Process %d: EigenVal: %f\n", my_rank, recvBuff);

    subtractEigenValue(m1, recvBuff);

    //MtxDisplay(m1);
    MtxToReducedREForm(m1);
    EL_Type* eigen_vector = (EL_Type*)malloc(sizeof(EL_Type) * (m1->dim_x - 1));
    EigenVectorNormalize(m1, eigen_vector, recvBuff);
    //MtxDisplay(m1);
    //printVec(eigen_vector, m1->dim_x - 1);
    
    EL_Type* gather_buff = NULL;
    EL_Type* u = NULL;
    if(my_rank == 0){
        gather_buff = (EL_Type*)TALLOC(MATRIX_SIZE * MATRIX_SIZE, EL_Type*);
	u = (EL_Type*)TALLOC(MATRIX_SIZE * MATRIX_SIZE, EL_Type*);
    }
    
    EL_Type *u_components = (EL_Type*)TALLOC(m1->dim_x-1, EL_Type*);
    MtxMulVec(m1, eigen_vector, u_components);

    
    MPI_Gather(eigen_vector, m1->dim_x-1, MPI_DOUBLE, gather_buff, m1->dim_x-1, \
	       MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(u_components, m1->dim_y, MPI_DOUBLE, u, m1->dim_y, MPI_DOUBLE, 0,\
	       MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    /*if(my_rank == 0) {
      for (i = 0; i < MATRIX_SIZE; i++) {
	printf("--------------------\nEigenvector %d\n----------------------\n", i);
	for (j = 0; j < MATRIX_SIZE; j++) {
	  
	  printf("%f\n", gather_buff[i+j]);
	}
      }
      }*/
    if(my_rank == 0) printf("parallel  time: %f\n", k, omp_get_wtime() - totalTime);
    return 0;
}

void computeQR(bool verbose, bool definition, float **matrixR, float **matrixQ, int rank, int lines)
{
    bool root = false;
    int i,j,size;
    int displs[MAX_PROCESSES], send_counts[MAX_PROCESSES], displs2[MAX_PROCESSES],\
      send_counts2[MAX_PROCESSES];
    float **mat,**p,**matTmp,**matTmp2;
    float *vec, coef;
    if(rank == 0) root = true;
    
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if(rank == 0) root = true;
    
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    
    for(int i = 0; i < lines; i++){
        matrixQ[i][i] = 1;
    }

    int k,l,m,tmpLines,tmpLines2;
    

    for(i = 0; i < lines; i++){
        tmpLines = (lines-i)/size;
        if(rank == (size-1) && size > 1) tmpLines += (lines-i)%size;
        tmpLines2 = (lines)/size;
        if(rank == (size-1) && size > 1) tmpLines2 += (lines)%size;
        //tmp matrix for parallel computing
        if(i > 0 && i < lines-size){
            delete [] mat[0];
        }
        mat = create_matrix(tmpLines,lines-i);

        MPI_Bcast((void *)&matrixR[0][0], lines*lines, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast((void *)&matrixQ[0][0], lines*lines, MPI_FLOAT, 0, MPI_COMM_WORLD);

        // nastaveni pro posilani matice o rozmeru [lines-i][lines-i]
        int piece = (lines-i)/size;
        int radek = lines-i;
        send_counts[0] = (piece)*(radek);
        displs[0] = 0;
        if(size > 1) {
            for(k = 1; k < size-1; k++){
                send_counts[k] = piece*(radek);
                displs[k] = displs[k-1] + piece*(radek);
            }
            displs[size-1] = displs[size-2] + piece*(radek);
            send_counts[size-1] = (((lines-i)*(radek)) - displs[size-1]);
        }

        // nastaveni pro posilani matice o rozmeru [lines][lines-i]
        piece = (lines)/size;
        radek = lines-i;
        send_counts2[0] = (piece)*(radek);
        displs2[0] = 0;
        if(size > 1){
            for(k = 1; k < size-1; k++){
                send_counts2[k] = piece*(radek);
                displs2[k] = displs2[k-1] + piece*(radek);
            }
            displs2[size-1] = displs2[size-2] + piece*(radek);
            send_counts2[size-1] = (((lines)*(radek)) - displs2[size-1]);
        }

        float x = 0;
        if(i > 0 && i < lines-size){
            delete [] vec;
            vec = NULL;
        }
        vec = new float[lines-i];
        if(rank == 0){
            for(j = 0; j < lines-i; j++){
                vec[j] = -matrixR[j+i][i];
                x += vec[j]*vec[j];
            }
            
            x = sqrt(x);
            
            
            if(vec[0] > 0) x = -x;
            vec[0] = vec[0] + x;
            x = 0;
            for(j = 0; j < lines-i; j++){
                x += vec[j]*vec[j];
            }
            x = sqrt(x);
        }
        MPI_Bcast((void *)&x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Bcast((void *)&vec[0], lines-i, MPI_FLOAT, 0, MPI_COMM_WORLD);
        if(x > 0){
            if(rank == 0){
                //normalizovat vec
                for(j = 0; j < lines-i; j++){
                    vec[j] /= x;
                }
            }
            MPI_Bcast((void *)&vec[0], lines-i, MPI_FLOAT, 0, MPI_COMM_WORLD);
            //sestavit matici P
            if(i > 0 && i < lines-size){
                delete [] p[0];
                p = NULL;
            }
            p = create_matrix(lines-i,lines-i);
            MPI_Scatterv(&p[0][0],send_counts,displs,MPI_FLOAT,&mat[0][0],tmpLines*(lines-i),MPI_FLOAT,0,MPI_COMM_WORLD);
            for(k = 0; k < send_counts[rank]/(lines-i); k++){
                for(l = 0; l < lines-i; l++){
                    if((k+(displs[rank]/(lines-i))) == l) mat[k][l] = 1 - 2*vec[k+displs[rank]/(lines-i)]*vec[l];
                    else mat[k][l] = -2*vec[k+displs[rank]/(lines-i)]*vec[l];
                    //if(mat[k][l] == 0) printf("%d -- \n",k+displs[rank]/(lines-i));
                }
            }
            MPI_Gatherv(&mat[0][0],send_counts[rank],MPI_FLOAT,&p[0][0],send_counts,displs,MPI_FLOAT,0,MPI_COMM_WORLD);
            MPI_Bcast((void *)&p[0][0], (lines-i)*(lines-i), MPI_FLOAT, 0, MPI_COMM_WORLD);
            //nasobeni matic (paralelizace)
            //R
            for(k = 0; k < send_counts[rank]/(lines-i); k++){
                for(l = 0; l < lines-i; l++){
                    float tm = 0;
                    for(m = i; m < lines; m++){
                        tm += p[k+displs[rank]/(lines-i)][m-i]*matrixR[m][l+i];
                    }
                    mat[k][l] = tm;
                }
            }
            if(i > 0 && i < lines-size){
                delete [] matTmp[0];
                matTmp = NULL;
            }
            matTmp = create_matrix(lines-i,lines-i);
            MPI_Gatherv(&mat[0][0],send_counts[rank],MPI_FLOAT,&matTmp[0][0],send_counts,displs,MPI_FLOAT,0,MPI_COMM_WORLD);
            if(rank == 0){
                for(k = i; k < lines; k++){
                    for(l = i; l < lines; l++){
                        matrixR[k][l] = matTmp[k-i][l-i];
                    }
                }
            }
            //Q
            if(i > 0 && i < lines-size){
                delete [] mat[0];
                mat = NULL;
            }
            mat = create_matrix(tmpLines2,lines-i);
            for(k = 0; k < send_counts2[rank]/(lines-i); k++){
                for(l = 0; l < lines-i; l++){
                    float tm = 0;
                    for(m = i; m < lines; m++){
                        tm += matrixQ[k+displs2[rank]/(lines-i)][m]*p[m-i][l];
                    }
                    mat[k][l] = tm;
                }
            }
            if(i > 0 && i < lines-size){
                delete [] matTmp[0];
                matTmp = NULL;
            }
            matTmp = create_matrix(lines,lines-i);
            MPI_Gatherv(&mat[0][0],send_counts2[rank],MPI_FLOAT,&matTmp[0][0],send_counts2,displs2,MPI_FLOAT,0,MPI_COMM_WORLD);
            if(rank == 0){
                for(k = 0; k < lines; k++){
                    for(l = i; l < lines; l++){
                        matrixQ[k][l] = matTmp[k][l-i];
                    }
                }
            }
        }
    }
    return;
}
