

#ifndef UPA_MYMATH_H
#define UPA_MYMATH_H

#include <stdexcept>

namespace upa {

    //! @brief Integer power
    int ipow(int base, int exp);

    //! @brief Determinant for small matrices
    double det(int dim, const double* A);

    //! @brief Determinant for 2x2 matrices
    double det2x2(const double* A);

    //! @brief Inverse for small matrices
    void inverse(int dim, const double* A, double* Ainv);

    //! @brief Inverse for 2x2
    void inverse2x2(const double* A, double* Ainv);





}

#endif //UPA_MYMATH_H
