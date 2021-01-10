

#ifndef UPA_DIRECTSOLVER_H
#define UPA_DIRECTSOLVER_H

#include "sparse_CSR.h"

namespace upa {

    class DirectSolver {

    protected:
        /// Inputs
        int n;
        Sparse_CSR *A;
        double *b;

        /// Parameters
        int verbose;

        /// Outputs
        double *sol;

    public:

        /// Constructor & Destructor
        DirectSolver() = default;
        virtual ~DirectSolver() = default;

        /// Setup
        void setVerbosity(int verbosity);

        /// Solve
        virtual void solve(double *x0) = 0;

        /// Getters
        virtual void getSolution(double *x);

    };

}

#endif //UPA_DIRECTSOLVER_H
