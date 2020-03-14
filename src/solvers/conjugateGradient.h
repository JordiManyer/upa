//
// Created by bscuser on 3/6/20.
//

#ifndef UPA_COJUGATEGRADIENT_H
#define UPA_COJUGATEGRADIENT_H


class conjugateGradient {

private:
    /// Inputs
    int n;
    double* A;
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

    /// Constructor
    conjugateGradient(int problemSize, double* problemMatrix, double* problemVector);

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
