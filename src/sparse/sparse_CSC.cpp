

#include "sparse_CSC.h"

sparse_CSC::sparse_CSC(int dimension, int nonzeros) {
    n = dimension;
    nnz = nonzeros;

    cols = new int[n+1];
    rows = new int[nnz];
    values = new double[nnz];
}

double sparse_CSC::get(int i, int j) {
    for (int k = cols[j]; k < cols[j+1]; ++k) {
        if ( rows[k] == i) return values[k];
    }
    return 0.0;
}

void sparse_CSC::matvec(double* x, double* y) {
    for (int i = 0; i < n; ++i) y[i] = 0.0;

    for (int i = 0; i < n; ++i) {
        for (int k = cols[i]; k < cols[i+1]; ++k) {
            y[rows[k]] += values[k] * x[i];
        }
    }
}

