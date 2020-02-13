

#ifndef LINEARSOLVERS_LINEARSOLVERS_H
#define LINEARSOLVERS_LINEARSOLVERS_H

#include <cmath>

// LU decomposition with partial pivoting
void LU(int n, double* A, int* P);
void solveLU(int n, double* LU , int* P , double* b, double* x);

// Cholesky decomposition
void Chol(int n, double* M, double* L);
void solveChol(int n, double* L, double* b, double* x);

#endif //LINEARSOLVERS_LINEARSOLVERS_H
