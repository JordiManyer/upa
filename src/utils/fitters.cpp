

#include <iostream>
#include "fitters.h"
#include "../solvers/directSolvers.h"

double phi(double x, double rho);
void Basis(double x , int degree, double* P);
void DBasis(double x , int degree, double* P);



/********************************************************************
 *******************     PUBLIC FUNCTIONS    ************************
 ********************************************************************/

/// Polynomial Least squares 1D
void PLS1D( data1D* data_in, int degree, double* c) {
    int m = degree + 1; // Number of coeffs
    double M[m*m], f[m]; // System matrices
    double P[m], L[m*m];

    // Pointers to data structures
    int n; double *x, *y;
    n = data_in->n; x = data_in->x; y = data_in->y;

    // Init system to zero
    for (int i = 0; i < m; ++i) {
        f[i] = 0.0;
        for (int j = 0; j < m; ++j) M[i*m + j] = 0.0;
    }

    // Assembly of the system
    for (int k = 0; k < n; ++k) {
        Basis(x[k], degree, P); // Fill polynomials
        for (int i = 0; i < m; ++i) {
            f[i] += y[k] * P[i];
            for (int j = 0; j < m; ++j) M[i*m + j] += P[i]*P[j];
        }
    }

    // Solve system (Symmetric, positive definite)
    Chol(m, M, L);
    solveChol(m, L, f, c); // c holds the coefficients
}

/// Polynomial Least squares 2D
void PLS2D() {

}

/// Moving least squares 1D
void MLS1D(data1D* data_in, int degree, double rho, data1D* data_out) {
    int iL, iR;
    int m = degree + 1; // Number of coeffs
    double M[m*m], f[m]; // System matrices
    double P[m], L[m*m], c[m];
    double phik, hk;

    // Pointers to data structures
    int n_in; double *x_in, *y_in;
    n_in = data_in->n; x_in = data_in->x; y_in = data_in->y;
    int n_out; double *x_out, *y_out;
    n_out = data_out->n; x_out = data_out->x; y_out = data_out->y;

    // Loop in output points
    iL = 0; iR = 1;
    for (int q = 0 ; q < n_out ; ++q) {
        // Define characteristic length
        if (q == 0) hk = fabs( x_out[q] - x_out[q+1]);
        else if (q == n_out - 1) hk = fabs( x_out[q-1] - x_out[q]);
        else hk = fabs( x_out[q-1] - x_out[q+1]) / 2.0;
        hk *= rho;

        // Find biggest iL and lowest iR such that
        // x_in[iL] <  x_out[q] - hk  <  x_out[q] + hk  <  x_in[iR]
        while(iL < n_in and x_in[iL] < x_out[q] - hk) ++iL;
        while(iR < n_in and x_in[iR] < x_out[q] + hk) ++iR;
        if (iL != 0) --iL;

        // Init system to zero
        for (int i = 0; i < m; ++i) {
            f[i] = 0.0;
            for (int j = 0; j < m; ++j) M[i*m + j] = 0.0;
        }

        // Assembly of the system considering neighbours
        for (int k = iL; k < iR+1; ++k) {
            phik = phi(x_in[k]-x_out[q],hk);
            Basis(x_in[k], degree, P); // Fill polynomials
            for (int i = 0; i < m; ++i) {
                f[i] += y_in[k] * P[i] * phik;
                for (int j = 0; j < m; ++j) M[i*m + j] += P[i] * P[j] * phik;
            }
        }

        // Solve system (Symmetric, positive definite)
        Chol(m, M, L);
        solveChol(m, L, f, c); // c holds the coefficients

        // Solution for this point
        Basis(x_out[q], degree, P);
        y_out[q] = 0.0;
        for (int k = 0; k < m; ++k) y_out[q] += P[k] * c[k];
    }
}




/********************************************************************
 *******************     PRIVATE FUNCTIONS    ***********************
 ********************************************************************/
double phi(double x, double rho) {
    double r = fabs(x/rho);
    if (r >= 1.0) return 0;
    else if (r >= 0.5) return (4.0/3.0) * (1-r) * (1-r) * (1-r);
    else return (2.0/3.0) - 4.0 * r * r * (1-r);
}

/// Polynomial basis functions and their derivatives
void Basis(double x , int degree, double* P) {
    if (degree == 0) {
        P[0] = 1;
    } else if (degree == 1) {
        P[0] = 1;
        P[1] = x;
    } else if (degree == 2) {
        P[0] = 1;
        P[1] = x;
        P[2] = x*x;
    } else {
        std::cerr << " ERROR from FITTERS: Degree not implemented!!" << std::endl;
    }
}
void DBasis(double x , int degree, double* P) {
    if (degree == 0) {
        P[0] = 0;
    } else if (degree == 1) {
        P[0] = 0;
        P[1] = 1;
    } else if (degree == 2) {
        P[0] = 0;
        P[1] = 1;
        P[2] = 2.0*x;
    } else {
        std::cerr << " ERROR from FITTERS: Degree not implemented!!" << std::endl;
    }
}



