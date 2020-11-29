

#ifndef UPA_SOLVER_H
#define UPA_SOLVER_H

#include "sparse_CSR.h"


namespace upa {

    class Solver {

    protected:
        /// Inputs
        int n;
        Sparse_CSR *A;
        double *b;

        /// Parameters
        int verbose;
        double tol;

        /// Outputs
        bool hasConverged;
        double *sol;
        double *residual;
        int numIter;
        double finalError;

    public:

        /// Destructor
        virtual ~Solver() = default;

        /// Setup
        virtual void configure(double tolerance, int beVerbose) = 0;

        /// Solve
        virtual void solve(double *x0) = 0;

        /// Getters
        virtual bool getConvergence() = 0;
        virtual void getSolution(double *x) = 0;
        virtual void getResidual(double *x) = 0;
        virtual int getNumIter() = 0;
        virtual double getError() = 0;

    };

}

#endif //UPA_SOLVER_H
