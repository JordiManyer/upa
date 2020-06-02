

#ifndef UPA_SPARSE_CSR_H
#define UPA_SPARSE_CSR_H

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


class sparse_CSR {

public:
    // VARIABLES
    int n;          // Size of the matrix
    int nnz;        // Number of non-zero elements
    int* rows;      // Size n, list of first indexes
    int* cols;      // Size nnz, indexes
    double* values; // Size nnz, values

    // METHODS
    sparse_CSR(int dimension, int nonzeros);

    double get(int i, int j);
    void matvec(double* x, double* y);
};


#endif //UPA_SPARSE_CSR_H
