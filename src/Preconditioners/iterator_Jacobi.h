

#ifndef UPA_ITERATOR_JACOBI_H
#define UPA_ITERATOR_JACOBI_H

#include "iterator.h"
#include "sparse_CSR.h"

namespace upa {

    /**  Jacobi iterator
     *   Each iteration performs the action  x <- (1-w) x + w[x + D^-1 (b - Ax)]
     *   where
     *      1) w is a relaxation parameter (tuneable).
     *      2) A is the matrix system, D is its diagonal and b is the system vector.
     */
    class Iterator_Jacobi: public Iterator {

    public:
        Iterator_Jacobi(Sparse_CSR* problemMatrix, double relaxationParameter);
        ~Iterator_Jacobi() override = default;
        void iterate(double* x_in, double* b, double* x_out) override;

    private:
        int n;
        Sparse_CSR* A;
        double* D;
        double w;

    };

}


#endif //UPA_ITERATOR_JACOBI_H
