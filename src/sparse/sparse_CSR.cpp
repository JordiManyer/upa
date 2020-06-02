

#include "sparse_CSR.h"

sparse_CSR::sparse_CSR(int dimension, int nonzeros) {
    n = dimension;
    nnz = nonzeros;

    rows = new int[n+1];
    cols = new int[nnz];
    values = new double[nnz];
}


double sparse_CSR::get(int i, int j) {
    for (int k = rows[i]; k < rows[i+1]; ++k) {
        if ( cols[k] == j) return values[k];
    }
    return 0.0;
}

void sparse_CSR::matvec(double* x, double* y) {
    for (int i = 0; i < n; ++i) {
        y[i] = 0.0;
        for (int k = rows[i]; k < rows[i+1]; ++k) {
            y[i] += values[k] * x[cols[k]];
        }
    }
}

