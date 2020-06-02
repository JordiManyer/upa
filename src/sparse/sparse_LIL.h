#ifndef UPA_SPARSE_LIL_H
#define UPA_SPARSE_LIL_H

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

#include <algorithm>
#include <vector>
#include <utility>


class sparse_LIL {

public:
    // VARIABLES
    int n;
    int nnz;
    int max_col_size;
    std::vector <std::pair <int,double>> *A;

    // METHODS
    sparse_LIL(int number_of_rows);

    void assemble_elem(int i, int j, double value);

    double get_elem_by_index(int i, int j);
    double get_elem_by_position(int i, int pos_j);

    int get_elem_index(int i, int pos_j);
};

#endif //UPA_SPARSE_LIL_H
