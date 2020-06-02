
#include "conjugateGradient.h"
#include <math.h>
#include <iostream>

/*****************************************************************************************
 *************************              CONSTRUCTORS             *************************
 *****************************************************************************************/

// Dense Constructor
conjugateGradient::conjugateGradient(int problemSize, double* problemMatrix, double* problemVector) {

    // Inputs
    n = problemSize;
    A = problemMatrix;
    b = problemVector;

    // Outputs
    hasConverged = false;
    sol = new double[n];
    residual = new double[n];
    numIter = 0;
    finalError = -1;

    // Set default parameters
    configure(1.e-12,false,0);
}

// Sparse Constructor
conjugateGradient::conjugateGradient(int problemSize, sparse_CSR* problemMatrix, double* problemVector) {

    // Inputs
    n = problemSize;
    AS = problemMatrix;
    b = problemVector;

    // Outputs
    hasConverged = false;
    sol = new double[n];
    residual = new double[n];
    numIter = 0;
    finalError = -1;

    // Set default parameters
    configure(1.e-12,true,0);
}


/*****************************************************************************************
 *************************        OPTIONS AND PARAMETERS         *************************
 *****************************************************************************************/
void conjugateGradient::configure(double tolerance, bool isSparse, int beVerbose){
    sparse  = isSparse;
    verbose = beVerbose;
    tol     = tolerance;
}


/*****************************************************************************************
 *************************       CONJUGATE GRADIENT METHOD       *************************
 *****************************************************************************************/

void conjugateGradient::solveDense(double* x0){
    int k;
    double rk[n], rkold[n], dk[n], Adk[n], xk[n];
    double alpha, beta, err, aux;

    /// Initialisation
    k = 0; err = 0.0;
    for (int i = 0; i < n; ++i) {
        Adk[i] = 0.0;
        for (int j = 0; j < n; ++j) Adk[i] += A[i*n + j]*x0[j];

        xk[i] = x0[i];
        rk[i] = b[i] - Adk[i];
        dk[i] = rk[i];
        err += rk[i]*rk[i];
    }
    err = sqrt(err);

    /// Main Loop :
    while ( err > tol and k < 2*n ) {
        // Calculate A * dk
        for (int i = 0; i < n; ++i) {
            Adk[i] = 0.0;
            for (int j = 0; j < n; ++j) Adk[i] += A[i * n + j] * dk[j];
        }

        // alpha = (dk·rk) / (dk·A·dk)
        alpha = 0.0; aux = 0;
        for (int i = 0; i < n; ++i) {
            alpha += dk[i]*rk[i];
            aux   += dk[i]*Adk[i];
        }
        alpha /= aux;

        // xkp1 = xk + alpha * dk ; rkp1 = rk - alpha * A * dk
        for (int i = 0; i < n; ++i) {
            xk[i] += alpha * dk[i];

            rkold[i] = rk[i];
            rk[i] -= alpha * Adk[i];
        }

        // beta = (rkp1·rkp1) / (rk·rk)
        beta = 0.0; aux = 0.0;
        for (int i = 0; i< n; ++i) {
            beta += rk[i]*rk[i];
            aux  += rkold[i]*rkold[i];
        }
        beta /= aux;

        // dkp1 = rkp1 + beta * dk
        for (int i = 0; i < n; ++i) dk[i] = rk[i] + beta * dk[i];

        // Update stop conditions
        err = 0.0; double orthogonal = 0;
        for (int i = 0; i < n; ++i) {
            orthogonal += rk[i]*rkold[i];
            err += rk[i] * rk[i];
        }
        ++k;
    }

    /// Save outputs:
    if (k == n+1) hasConverged = false;
    else hasConverged = true;

    for (int i = 0; i < n; ++i) {
        sol[i] = xk[i];
        residual[i] = rk[i];
    }
    numIter = k;
    finalError = err;
}


void conjugateGradient::solveSparse(double* x0) {
    int k;
    double rk[n], rkold[n], dk[n], Adk[n], xk[n];
    double alpha, beta, err, aux;

    /// Initialisation
    k = 0; err = 0.0;
    AS->matvec(x0, Adk);
    for (int i = 0; i < n; ++i) {
        xk[i] = x0[i];
        rk[i] = b[i] - Adk[i];
        dk[i] = rk[i];
        err += rk[i]*rk[i];
    }
    err = sqrt(err);

    /// Main Loop :
    while ( err > tol and k < 2*n ) {
        // Calculate A * dk
        AS->matvec(dk, Adk);

        // alpha = (dk·rk) / (dk·A·dk)
        alpha = 0.0; aux = 0;
        for (int i = 0; i < n; ++i) {
            alpha += dk[i]*rk[i];
            aux   += dk[i]*Adk[i];
        }
        alpha /= aux;

        // xkp1 = xk + alpha * dk ; rkp1 = rk - alpha * A * dk
        for (int i = 0; i < n; ++i) {
            xk[i] += alpha * dk[i];

            rkold[i] = rk[i];
            rk[i] -= alpha * Adk[i];
        }

        // beta = (rkp1·rkp1) / (rk·rk)
        beta = 0.0; aux = 0.0;
        for (int i = 0; i< n; ++i) {
            beta += rk[i]*rk[i];
            aux  += rkold[i]*rkold[i];
        }
        beta /= aux;

        // dkp1 = rkp1 + beta * dk
        for (int i = 0; i < n; ++i) dk[i] = rk[i] + beta * dk[i];

        // Update stop conditions
        err = 0.0; double orthogonal = 0;
        for (int i = 0; i < n; ++i) {
            orthogonal += rk[i]*rkold[i];
            err += rk[i] * rk[i];
        }
        ++k;
    }

    /// Save outputs:
    if (k == n+1) hasConverged = false;
    else hasConverged = true;

    for (int i = 0; i < n; ++i) {
        sol[i] = xk[i];
        residual[i] = rk[i];
    }
    numIter = k;
    finalError = err;
}


void conjugateGradient::solve(double* x0) {
    if (sparse) solveSparse(x0);
    else solveDense(x0);
}


/*****************************************************************************************
 *************************                GETTERS                *************************
 *****************************************************************************************/
bool   conjugateGradient::getConvergence() {
    return hasConverged;
}
void   conjugateGradient::getSolution(double* x) {
    for (int i = 0; i < n; ++i) x[i] = sol[i];
}
void   conjugateGradient::getResidual(double* r) {
    for (int i = 0; i < n; ++i) r[i] = residual[i];
}
int    conjugateGradient::getNumIter() {
    return numIter;
}
double conjugateGradient::getError() {
    return finalError;
}

