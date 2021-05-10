

#ifndef UPA_MYMATH_H
#define UPA_MYMATH_H

namespace upa {

    ///                             MISCELLANEOUS

    //! @brief Integer power
    int ipow(int base, int exp);


    ///                             MATRIX AND VECTOR MANIPULATION

    //! @brief Determinant for small matrices
    double det(int dim, const double* A);

    //! @brief Determinant for 2x2 matrices
    double det2x2(const double* A);

    //! @brief Inverse for small matrices
    void inverse(int dim, const double* A, double* Ainv);

    //! @brief Inverse for 2x2
    void inverse2x2(const double* A, double* Ainv);

    //! @brief Dot product
    double dot(int n, const double* x, const double* y);

    ///                              INTEGRATION

    //! @brief Performs the integral int{ y(x) dx} using Simpson's rule on the provided samples.
    double Simpson(int n, const double *x, const double *y);

    //! @brief Performs the integral int{ y(x) dx} using Simpson's rule on the samples iStart <= i < iEnd.
    double Simpson(int iStart, int iEnd, const double *x, const double *y);

    //! @brief Returns Gauss-Legendre quadratures in the interval [a,b]
    void GLquad(int n, double a, double b, double* w, double* z);


}

#endif //UPA_MYMATH_H
