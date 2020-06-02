

#include "sparse.h"

// Converts a LIL matrix to CSR
sparse_CSR* LILtoCSR(sparse_LIL* lil) {
    int k;
    sparse_CSR* csr = new sparse_CSR(lil->n,lil->nnz);
    for (int i = 0 ; i < lil->n ; ++i) {
        csr->rows[i+1] = csr->rows[i] + lil->A[i].size();

        k = csr->rows[i];
        for (auto it = lil->A[i].begin() ; it < lil->A[i].end() ; ++it) {
            csr->cols[k]   = it->first;
            csr->values[k] = it->second;
            ++k;
        }
    }
    return csr;
}