

#ifndef UPA_SOLVER_CG_H
#define UPA_SOLVER_CG_H

#include "sparse_CSR.h"
#include "solver.h"

namespace upa {

    class Solver_CG: public Solver {

    private:

    public:

        /// Constructors & Destructors
        Solver_CG(int problemSize, Sparse_CSR* problemMatrix, double *problemVector);
        ~Solver_CG() override = default;

        /// Modify default parameters
        void configure(double tolerance, int beVerbose) override;

        /// Solve the system using the conjugate gradient method
        void solve(double *x0) override;

        /// Getters
        bool getConvergence() override;
        void getSolution(double *x) override;
        void getResidual(double *x) override;
        int getNumIter() override;
        double getError() override;

    };

}

#endif //UPA_SOLVER_CG_H
