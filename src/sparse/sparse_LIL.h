

#ifndef UPA_SPARSE_LIL_H
#define UPA_SPARSE_LIL_H


#include <algorithm>
#include <vector>
#include <utility>

class sparse_LIL {

private:
    int n;
    int num_non_zeros;
    int max_col_size;
    std::vector <std::pair <int,double>> *A;

public:
    sparse_LIL(int number_of_rows);

    void assemble_elem(int i, int j, double value);

    double get_elem_by_index(int i, int j);
    double get_elem_by_position(int i, int pos_j);

    int get_elem_index(int i, int pos_j);
};


#endif //UPA_SPARSE_LIL_H
