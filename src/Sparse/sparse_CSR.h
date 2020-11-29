

#ifndef UPA_SPARSE_CSR_H
#define UPA_SPARSE_CSR_H

#include <vector>
#include "sparse_LIL.h"

namespace upa {
    /*************************************************************************
     ****       Compressed Storage Row-major (CSR) matrix storage      *******
     *************************************************************************
     * The matrix is represented as three arrays:
     *  - rows : Array containing at its i-th position the start of the i-th
     *           row nonzero elements. I.e the nonzero elements of the i-th
     *           row are located at positions rows[i], ..., rows[i+1].
     *  - cols : Array containing the column indexes of the nonzero entries.
     *  - values : Array containing the values of the nonzero entries.
     */
    class Sparse_CSR {

    public:
        // VARIABLES
        int n;          // Size of the matrix
        int nnz;        // Number of non-zero elements
        int *rows;      // Size n, list of first indexes
        int *cols;      // Size nnz, indexes
        double *values; // Size nnz, values
        int *diag;      // Size n, list of indexes for the diagonal elements

        // METHODS
        Sparse_CSR(int matrix_size, int nonzeros);     // Constructor
        explicit Sparse_CSR(Sparse_LIL *lil);          // Constructor from LIL
        explicit Sparse_CSR(Sparse_LIL *lil,bool *select);
        ~Sparse_CSR() = default;

        double operator() (int row, int col);          // Access operator

        void matvec(double *x, double *y);
        void getDiag();
    };
}

#endif //UPA_SPARSE_CSR_H
