

#include "myMath.h"
#include <stdexcept>
#include <cmath>

const double PI = 3.14159265358979323846;

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


    void GLquad(int n, double a, double b, double *w, double *z) {
        n = n-1; int n1 = n+1; int n2 = n+2;
        double y[n1], y_old[n1], x, maxErr;
        double Pk[n1], Pkm1[n1], dPk[n1];

        if (n == 0) { // Special case
            w[0] = 2.0; z[0] = 0.0; return;
        }

        for (int i = 0; i < n1; ++i) {
            y_old[i] = 2.0;
            x = -1.0 + 2.0*double(i)/double(n); // Equally spaced points in [-1,1]
            y[i] = cos((2.0*i+1.0)*PI/(2.0*n+2.0)) + (0.27/n1) * sin(PI*x*n/n2); // Initial guess
        }

        maxErr = 10.0;
        while ( maxErr > 1.e-10 ) { // Newton-Raphson loop
            for (int i = 0; i < n1; ++i) { Pk[i] = y[i]; Pkm1[i] = 1.0;} // Setup P0 = 1.0, P1 = y
            // Compute P_{n1} by using the Legendre recurrence   k·P_{k} = (2k-1)·y·P_{k-1} + (k-1)·P_{k-2}
            for (int k = 2; k < n2; ++k) {
                for (int i = 0; i < n1; ++i) {
                    x       = Pk[i]; // Save value
                    Pk[i]   = ((2.0*k-1.0)*y[i]*Pk[i] - (k-1.0)*Pkm1[i])/k; // Compute P_{k}
                    Pkm1[i] = x;
                }
            }
            // Compute the derivative of Pk
            for (int i = 0; i < n1; ++i) dPk[i] = n2 * (Pkm1[i]-y[i]*Pk[i]) / (1.0-y[i]*y[i]);
            // Compute next iterate and update convergence conditions
            maxErr = 0.0;
            for (int i = 0; i < n1; ++i) {
                y_old[i] = y[i];
                y[i] = y_old[i] - Pk[i]/dPk[i];
                if(fabs(y_old[i]-y[i]) > maxErr) maxErr = fabs(y_old[i]-y[i]);
            }
        }

        // Create quadratures
        for (int i = 0; i < n1; ++i) {
            //  Linear map from[-1,1] to [a,b]
            z[i] = (a*(1.0-y[i])+b*(1.0+y[i]))/2.0;
            //  Compute the weights
            w[i] = (b-a) / ((1.0-y[i]*y[i]) * dPk[i]*dPk[i]) * double(n2*n2)/double(n1*n1);
        }
    }




}