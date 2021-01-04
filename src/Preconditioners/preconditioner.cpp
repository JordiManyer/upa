
#include "preconditioner.h"


namespace upa {

    // Constructor for a general preconditioner
    preconditioner::preconditioner() {
    }

    // Build the Jacobi preconditioner
    void preconditioner::jacobi(Sparse_CSR *A) {
        type = ""
        "jacobi";
        n = A->n;

        // For Jacobi, we have M = diag(A)
        D = new double[n];
        bool found;
        for (int i = 0; i < n; ++i) {
            found = false;
            for (int k = A->rows[i]; not found and k < A->rows[i + 1]; ++k) {
                if (A->cols[k] == i) {
                    D[i] = A->values[k];
                    found = true;
                }
            }
            if (found == false) D[i] = 0; // WARNING : Output error!!
        }
    }

    // Build the SSOR preconditioner
    void preconditioner::ssor(Sparse_CSR *A, double w) {
        type = "ssor";
        n = A->n;

        // M = LU,   with   L = I - w E inv(D)   and   U = D - w F
        M = new Sparse_CSR(A->n, A->nnz);
        bool found;
        for (int i = 0; i < n; ++i) {
            M->rows[i] = A->rows[i];

            // Find ith diagonal element (inverted)
            found = false;
            for (int k = A->rows[i]; not found and k < A->rows[i + 1]; ++k)
                if (A->cols[k] == i) {
                    M->diag[i] = k;
                    M->cols[k] = i;
                    M->values[k] = 1.0 / A->values[k];
                    found = true;
                }
            if (found == false) M->diag[i] = -1; // WARNING : Output error!!

            // Fill ith row
            for (int k = A->rows[i]; k < M->diag[i]; ++k) { // L = I - w E inv(D)
                M->cols[k] = A->cols[k];
                M->values[k] = 1.0 + w * A->values[k] * M->values[M->diag[i]];
            }
            for (int k = M->diag[i] + 1; k < A->rows[i + 1]; ++k) { // U = D - w F
                M->cols[k] = A->cols[k];
                M->values[k] = 1.0 / M->values[M->diag[i]] + w * A->values[k];
            }
        }
        M->rows[n] = A->rows[n];
    }


/***********************************************************************************
 *************************           SOLVERS          ******************************
 ***********************************************************************************/

    // General solver
    void preconditioner::solve(double *x, double *y) {
        if      (type == "jacobi") solve_jacobi (double *x, double *y);
        else if (type == "ssor")   solve_ssor   (double *x, double *y);
    }

    // Jacobi Solver
    void preconditioner::solve_jacobi(double *x, double *y) {

    }

    // SSOR Solver
    void preconditioner::solve_ssor(double *x, double *y) {

    }

}