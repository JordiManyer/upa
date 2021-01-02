
#ifndef UPA_DEBUGIO_H
#define UPA_DEBUGIO_H

#include <iostream>
#include <iomanip>

namespace upa {

    template <typename T>
    void printArray(size_t n, const T array[]) {
        std::cout << std::endl;
        for (size_t i = 0; i < n; ++i) std::cout << std::setw(6) << array[i] << " ";
        std::cout << std::endl;
    }

    template<typename T>
    void printMatrix(size_t n, size_t m, const T matrix[]) {
        std::cout << std::endl;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) std::cout << std::setw(6) << matrix[i*n+j] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

}


#endif //UPA_DEBUGIO_H
