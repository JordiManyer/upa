

#ifndef UPA_DENSESOLVERS_H
#define UPA_DENSESOLVERS_H

#include <cmath>

// LU decomposition with partial pivoting
void LU(int n, double* A, int* P);
void solveLU(int n, double* LU , int* P , double* b, double* x);

// Cholesky decomposition
void Chol(int n, double* M, double* L);
void solveChol(int n, double* L, double* b, double* x);

#endif //UPA_DENSESOLVERS_H
