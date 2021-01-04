

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
        int maxIter;

        /// Outputs
        bool hasConverged;
        double *sol;
        double *residual;
        int numIter;
        double finalError;

    public:

        /// Constructor & Destructor
        Solver() = default;
        virtual ~Solver() = default;

        /// Setup
        void setTolerance(double tolerance);
        void setIterations(int maxNumIter);
        void setVerbosity(int verbosity);

        /// Solve
        virtual void solve(double *x0) = 0;

        /// Getters
        virtual bool getConvergence();
        virtual void getSolution(double *x);
        virtual void getResidual(double *x);
        virtual int getNumIter();
        virtual double getError();

    };

}

#endif //UPA_SOLVER_H
