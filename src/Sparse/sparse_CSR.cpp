

#include "sparse_CSR.h"
#include <iostream>

namespace upa {


    // Constructor for a sparse CSR matrix
    Sparse_CSR::Sparse_CSR(int matrix_size, int nonzeros) {
        n = matrix_size;
        nnz = nonzeros;

        rows = new int[n + 1];
        cols = new int[nnz];
        values = new double[nnz];
    }

    // Constructs a CSR matrix from a LIL matrix
    Sparse_CSR::Sparse_CSR(Sparse_LIL *lil) {
        n = lil->n;
        nnz = lil->nnz;

        rows = new int[n + 1];
        cols = new int[nnz];
        values = new double[nnz];

        rows[0] = 0;
        int k = 0;
        for (int i = 0; i < lil->n; ++i) {
            rows[i + 1] = rows[i] + lil->A[i].size();

            k = rows[i];
            for (auto it = lil->A[i].begin(); it < lil->A[i].end(); ++it) {
                cols[k] = it->first;
                values[k] = it->second;
                ++k;
            }
        }
    }

    // Constructs a CSR matrix from a LIL matrix, selecting only certain rows and cols
    Sparse_CSR::Sparse_CSR(Sparse_LIL *lil, bool *select) {

        // Count selected rows and number of non-zero entries
        n = 0;
        nnz = 0;
        for (int i = 0; i < lil->n; ++i) {
            if (select[i]) { // Only consider selected rows
                n++;
                for (auto e = lil->A[i].begin(); e < lil->A[i].end(); ++e)
                    if (select[e->first]) ++nnz; // Only consider selected columns
            }
        }

        rows = new int[n + 1];
        cols = new int[nnz];
        values = new double[nnz];

        int iCsr = 0; int k = 0; int skipped = 0;
        for (int iLil = 0; iLil < lil->n; ++iLil) {
            if (select[iLil]) {
                rows[iCsr + 1] = rows[iCsr] + lil->A[iLil].size();

                // TODO: We need to modify the column... -> Count skipped dofs and substract from .first
                k = rows[iCsr];
                for (auto e = lil->A[iLil].begin(); e < lil->A[iLil].end(); ++e)
                    if (select[e->first]) {
                        cols[k] = e->first - skipped;
                        values[k] = e->second;
                        ++k;
                    }
                iCsr++;
            } else ++skipped;
        }
    }

    // Define operator (i,j) for matrix access
    double Sparse_CSR:: operator() (int i, int j) {
        for (int k = rows[i]; k < rows[i + 1]; ++k) {
            if (cols[k] == j) return values[k];
        }
        return 0.0;
    }

    // Left-Multiply matrix by x, result in y
    void Sparse_CSR::matvec(double *x, double *y) {
        for (int i = 0; i < n; ++i) {
            y[i] = 0.0;
            for (int k = rows[i]; k < rows[i + 1]; ++k) {
                y[i] += values[k] * x[cols[k]];
            }
        }
    }

    // Fills a vector with the indexes of the diagonal elements
    void Sparse_CSR::getDiag() {
        diag = new int[n];
        bool found;
        for (int i = 0; i < n; ++i) {
            found = false;
            for (int k = rows[i]; k < rows[i + 1] and not found; ++k) {
                if (cols[k] == i) {
                    diag[i] = k;
                    found = true;
                }
            }
            if (not found) diag[i] = -1;
        }
    }

    void Sparse_CSR::getDiag(double* D) {
        bool found;
        for (int i = 0; i < n; ++i) {
            found = false;
            for (int k = rows[i]; k < rows[i + 1] and not found; ++k) {
                if (cols[k] == i) {
                    D[i] = values[k];
                    found = true;
                }
            }
            if (not found) D[i] = 0.0;
        }
    }

}

