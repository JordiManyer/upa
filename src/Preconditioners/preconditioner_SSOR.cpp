

#include "preconditioner_SSOR.h"


namespace upa {

    Preconditioner_SSOR::Preconditioner_SSOR(Sparse_CSR *A, double w) {
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
            if (not found) M->diag[i] = -1; // WARNING : Output error!!

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


    void Preconditioner_SSOR::apply(double *r_in, double *r_out) {

    }

}