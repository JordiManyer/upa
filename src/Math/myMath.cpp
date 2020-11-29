

#include "myMath.h"

namespace upa {

    int ipow(int base, int exp) {
        int res = 1;
        while (exp-- > 0) res *= base;
        return res;
    }


    double det(int dim, const double* A) {
        switch (dim) {
            case 1:
                return *A;
            case 2:
                return det2x2(A);
            default:
                throw std::runtime_error("upa_Math: Determinant not implemented.");
                return 0;

        }
    }

    double det2x2(const double* A) {
        return A[0]*A[3] - A[1]*A[2];
    }


    void inverse(int dim, const double* A, double* Ainv) {
        switch (dim) {
            case 1:
                *Ainv =  1.0/(*A); return;
            case 2:
                inverse2x2(A,Ainv); return;
            default:
                throw std::runtime_error("upa_Math: Inverse not implemented.");
                return;

        }
    }

    void inverse2x2(const double* A, double* Ainv) {
        double d = det2x2(A);
        Ainv[0] = A[3]/d;
        Ainv[1] = -A[1]/d;
        Ainv[2] = -A[2]/d;
        Ainv[3] = A[0]/d;
    }


}