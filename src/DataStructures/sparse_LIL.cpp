

#include "sparse_LIL.h"
#include <algorithm>
#include <utility>

namespace upa {

    Sparse_LIL::Sparse_LIL(int matrix_size) {
        n = matrix_size;
        nnz = 0;
        max_col_size = 0;
        A = new std::vector <std::pair<int, double>>[n];
    }

    double Sparse_LIL::operator() (int i, int j) {
        auto it = find_if(A[i].begin(), A[i].end(),
                     [&j](const std::pair<int, double> &element) { return element.first == j; });
        if (it == A[i].end()) return 0.0; // Element not in matrix (it is a zero)
        else return it->second;           // Element in matrix
    }

    void Sparse_LIL::put(int i, int j, double value) {
        auto it = find_if(A[i].begin(), A[i].end(),
                     [&j](const std::pair<int, double> &element) { return element.first == j; });

        if (it == A[i].end()) { // Element does not exist --> Create element
            A[i].push_back(std::make_pair(j, value));
            sort(A[i].begin(), A[i].end());
            ++nnz;
        } else { // Element already exists --> Add value to element
            it->second = value;
        }
    }

    void Sparse_LIL::assemble(int i, int j, double value) {
        auto it = find_if(A[i].begin(), A[i].end(),
                          [&j](const std::pair<int, double> &element) { return element.first == j; });

        if (it == A[i].end()) { // Element does not exist --> Create element
            A[i].push_back(std::make_pair(j, value));
            sort(A[i].begin(), A[i].end());
            ++nnz;
        } else { // Element already exists --> Add value to element
            it->second += value;
        }
    }

    void Sparse_LIL::removeRow(int i) {
        nnz -= A[i].size();
        A[i].erase(A[i].begin(),A[i].end());
    }

}

