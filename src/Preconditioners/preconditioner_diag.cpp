

#include "preconditioner_diag.h"

namespace upa {

    Preconditioner_Diag::Preconditioner_Diag(Sparse_CSR *A) {
        n = A->n;
        D = new double[n];
        A->getDiag(D);
    }


    void Preconditioner_Diag::apply(double *r_in, double *r_out) {
        for(int i = 0; i < n; ++i) r_out[i] = r_in[i]/D[i];
    }

}