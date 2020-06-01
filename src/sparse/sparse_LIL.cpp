

#include "sparse_LIL.h"


sparse_LIL::sparse_LIL(int number_of_rows) {
    n = number_of_rows;
    num_non_zeros = 0;
    max_col_size = 0;
    A = new std::vector <std::pair <int,double>> [n];
}

void sparse_LIL::assemble_elem(int i, int j, double value){
    std::vector<std::pair<int,double>>::iterator it = find_if( A[i].begin(), A[i].end(),
                                                              [&j](const std::pair<int, double>& element){ return element.first == j;} );

    if (it == A[i].end()) { // Element does not exist --> Create element
        A[i].push_back( std::make_pair(j,value) );
        sort(A[i].rbegin(), A[i].rend());
    } else { // Element already exists --> Add value to element
        it->second += value;
    }
}

double sparse_LIL::get_elem_by_index(int i, int j) {
    std::vector<std::pair<int,double>>::iterator it = find_if( A[i].begin(), A[i].end(),
                                                               [&j](const std::pair<int, double>& element){ return element.first == j;} );
    if (it == A[i].end()) return 0.0; // Element not in matrix (it is a zero)
    else return it->second;       // Element in matrix
}

double sparse_LIL::get_elem_by_position(int i, int pos_j) {
    return A[i][pos_j].second;
}

int sparse_LIL::get_elem_index(int i, int pos_j) {
    return A[i][pos_j].first;
}

