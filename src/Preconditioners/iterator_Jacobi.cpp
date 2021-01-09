

#include "iterator_Jacobi.h"


namespace upa {

    Iterator_Jacobi::Iterator_Jacobi(Sparse_CSR* problemMatrix, double relaxationParameter) {
        n = problemMatrix->n;
        A = problemMatrix;
        w = relaxationParameter;
        D = new double[n];
        A->getDiag(D);
    }

    void Iterator_Jacobi::iterate(double* x_in, double* b, double* x_out) {
        A->matvec(x_in, x_out);
        for (int i = 0; i < n; ++i) {
            x_out[i] = (b[i] - x_out[i])/D[i] + x_in[i];
            x_out[i] = w * x_out[i] + (1.0 - w) * x_in[i];
        }
    }

}