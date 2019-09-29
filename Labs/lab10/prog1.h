void serMM2CSR(int* I, int* J, double* valMM, int M, int N, int nz, double* valCSR, int* col_ind, int* row_ptr);

void serSpMatVec(double* val, int* col_ind, int* row_ptr, int M, int N, int nz, double* x, double* y);
