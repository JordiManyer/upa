

#ifndef UPA_SPARSE_CSC_H
#define UPA_SPARSE_CSC_H


namespace upa {

    // Compressed Storage Column-major
    class Sparse_CSC {

    public:
        // VARIABLES
        int n;          // Size of the matrix
        int nnz;        // Number of non-zero elements
        int* cols;      // Size n, list of first indexes
        int* rows;      // Size nnz, indexes
        double* values; // Size nnz, values

        // METHODS
        Sparse_CSC(int dimension, int nonzeros);

        double get(int i, int j);
        void matvec(double* x, double* y);

    };

}

#endif //UPA_SPARSE_CSC_