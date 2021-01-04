
#include "solver.h"

namespace upa {


    /*****************************************************************************************
    *************************        OPTIONS AND PARAMETERS         *************************
    *****************************************************************************************/
    void Solver::setTolerance(double tolerance) {
        tol = tolerance;
    }

    void Solver::setIterations(int maxNumIter) {
        maxIter = maxNumIter;
    }

    void Solver::setVerbosity(int verbosity) {
        verbose = verbosity;
    }



    /***********************************************************
    *********                GETTERS                ************
    ************************************************************/
    bool Solver::getConvergence() {
        return hasConverged;
    }

    void Solver::getSolution(double *x) {
        for (int i = 0; i < n; ++i) x[i] = sol[i];
    }

    void Solver::getResidual(double *r) {
        for (int i = 0; i < n; ++i) r[i] = residual[i];
    }

    int Solver::getNumIter() {
        return numIter;
    }

    double Solver::getError() {
        return finalError;
    }


}
