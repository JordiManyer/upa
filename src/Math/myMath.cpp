

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


    double dot(int n, const double *x, const double *y) {
        double res = 0.0;
        for (int i = 0; i < n; ++i) res += x[i]*y[i];
        return res;
    }


    double Simpson(int n, const double *x, const double *y) {
        double res = 0.0; int i = 0;
        while (i < n-2) {
            res += (y[i] + 4 * y[i+1] + y[i+2]) * (x[i+2]-x[i]) / 6.0;
            i = i+2;
        }
        // If n != 2*k + 1, use trapezoidal rule on last segment
        if (i != n) res += (y[n-2] + y[n-1]) * (x[n-1]-x[n-2]) / 2.0;
        return res;
    }


    double Simpson(int iStart, int iEnd, const double *x, const double *y) {
        double res = 0.0; int i = iStart;
        while (i < iEnd-2) {
            res += (y[i] + 4 * y[i+1] + y[i+2]) * (x[i+2]-x[i]) / 6.0;
            i = i+2;
        }
        // If n != 2*k + 1, use trapezoidal rule on last segment
        if (i != iEnd) res += (y[iEnd-2] + y[iEnd-1]) * (x[iEnd-1]-x[iEnd-2]) / 2.0;
        return res;
    }


}