

#include <iostream>
#include "linearSolvers.h"

/***********************************************************************
 **************  LU DECOMPOSITION WITH PARTIAL PIVOTING  ***************
 ***********************************************************************/

/// LU decomposition with row pivoting: P·A = L·U
///    · LU is stored into A, P is given as a permutation vector
///    · Works for any non-singular matrix. Not precise for near-singular cases.
///    · WARNING: A is destroyed in the process.
void LU(int n, double* A, int* P) {
    int pivot;
    double sum, max, aux;

    // For all columns
    for (int k = 0; k < n ; ++k) {
        // Choose pivot
        pivot = k; max = fabs(A[k*n + k]);
        for (int j = k; j < n; ++j) if (max < fabs(A[k*n + j])) {
            max = fabs(A[k*n + j]);
            pivot = j;
        }
        P[k] = pivot;

        // If needed, issue conditioning warning:
        if (max < 1.e-3) std::cout << "WARNING from LU -> Matrix is bad conditioned!" << std::endl;

        // If needed, swap rows
        if (pivot != k) for (int j = 0; j < n; ++j) {
                aux = A[k*n + j];
                A[k*n + j] = A[pivot*n + j];
                A[pivot*n + j] = aux;
        }

        // We calculate the kth row of U
        // Ukj = Akj - sum_{i=1,...k-1} Uij Lki  , j >= k
        for (int j = k; j < n; ++j) {
            sum = 0.0;
            for (int i = 0; i < k; ++i) sum += A[i*n + j]*A[k*n + i];
            A[k*n + j] -= sum;
        }

        // We calculate the kth column of L
        // Lik = ( Aik - sum_{j=1,...k-1} Ujk Lij ) / Ujj , k < i
        // Lkk = 1.0 implicitly
        for (int i = k+1; i < n; ++i) {
            sum = 0.0;
            for (int j = 0; j < k; ++j) sum += A[j*n + k]*A[i*n + j];
            A[i*n + k] -= sum;
            A[i*n + k] /= A[k*n + k];
        }
    }

}


/// Solves A·x = b, with P·A = L·U it's LU decomposition with partial pivoting.
/// Pivoting: A·x = b -> P·A·x = P·b -> L·U·x = P·b is the system to be solved
void solveLU(int n, double* LU , int* P , double* b, double* x) {
    double sum, y[n];

    // First, we solve L·y = P·b (Lower triangular system, diagonal assumed to be 1.0)
    for (int i = 0; i < n; ++i) {
        sum = 0.0;
        for (int j = 0; j < i ; ++j) sum += LU[i*n+j] * y[j];
        y[i] = b[P[i]] - sum ;
    }

    // Finally, we solve U·x = y (Upper triangular system)
    for (int i = n-1; i >= 0; --i) {
        sum = 0.0;
        for (int j = i+1; j < n ; ++j) sum += LU[i*n+j] * x[j];
        x[i] = (y[i] - sum) / LU[i*n+i];
    }
}



/***********************************************************************
 ******************  CHOLESKY DECOMPOSITION  ***************************
 ***********************************************************************/

/// Decomposes M = L · L', with L a lower triangular matrix.
///    · M is not destroyed in the process
///    · WARNING: Cholesky only works for symmetric positive-definite matrixes.
void Chol(int n, double* M, double* L) {
    double sum;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            sum = 0.0;
            if (j == i) { // Diagonal elements
                for (int k = 0; k < j; ++k)
                    sum += pow(L[n*j+k], 2);
                L[n*j+j] = sqrt(M[n*j+j] - sum);
            } else { // Non-diagonal elements
                for (int k = 0; k < j; ++k)
                    sum += (L[n*i+k] * L[n*j+k]);
                L[n*i+j] = (M[n*i+j] - sum) / L[n*j+j];
            }
        }
    }
}

/// Solves M·x = L·L'·x = b by using the Cholesky decomposition of M = L·L'
void solveChol(int n, double* L, double* b, double* x) {
    double sum, y[n];

    // First, we solve L·y = b (Lower triangular system)
    for (int i = 0; i < n; ++i) {
        sum = 0.0;
        for (int j = 0; j < i ; ++j) sum += L[i*n+j] * y[j];
        y[i] = (b[i] - sum) / L[i*n+i];
    }

    // Finally, we solve L'·x = y (Upper triangular system)
    for (int i = n-1; i >= 0; --i) {
        sum = 0.0;
        for (int j = i+1; j < n ; ++j) sum += L[j*n+i] * x[j];
        x[i] = (y[i] - sum) / L[i*n+i];
    }
}

