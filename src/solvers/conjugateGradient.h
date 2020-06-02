

#ifndef UPA_COJUGATEGRADIENT_H
#define UPA_COJUGATEGRADIENT_H

#include "sparse/sparse.h"

class conjugateGradient {

private:
    /// Inputs
    int n;
    double* A;
    sparse_CSR* AS;
    double* b;

    /// Parameters
    bool sparse;
    int verbose;
    double tol;

    /// Outputs
    bool    hasConverged;
    double* sol;
    double* residual;
    double  numIter;
    double  finalError;

    /// Private functions
    void solveDense(double* x0);
    void solveSparse(double* x0);

public:

    /// Constructors
    conjugateGradient(int problemSize, double* problemMatrix, double* problemVector); // dense
    conjugateGradient(int problemSize, sparse_CSR* problemMatrix, double* problemVector); //sparse

    /// Modify default parameters
    void configure(double tolerance, bool isSparse = false, int beVerbose = 0);

    /// Solve the system using the conjugate gradient method
    void solve(double* x0);

    /// Getters
    bool   getConvergence();
    void   getSolution(double* x);
    void   getResidual(double* x);
    int    getNumIter();
    double getError();

};


#endif //UPA_COJUGATEGRADIENT_H
