#ifndef UPA_SPARSE_LIL_H
#define UPA_SPARSE_LIL_H

#include <vector>

namespace upa {

    /*************************************************************************
     **********         List of Lists (LIL) matrix storage         ***********
     *************************************************************************
     * The matrix is represented as an array of std::vector of pairs in the
     * following way :
     * Each row 'i' has a vector A[i] of dynamics dimension which holds the
     * nonzero values present on the i-th row.
     * Elements of the vector A[i] are pairs<int,double> where the integer
     * contains the column of the entry and the double contains the matrix entry.
     */
    class Sparse_LIL {

    public:
        // VARIABLES
        int n;
        int nnz;
        int max_col_size;
        std::vector <std::pair<int, double>> *A;

        // METHODS
        Sparse_LIL(int matrix_size);
        ~Sparse_LIL() = default;

        double operator() (int i, int j);

        void put(int i, int j, double value);
        void assemble(int i, int j, double value);

        void removeRow(int i);

    };
}
#endif //UPA_SPARSE_LIL_H
