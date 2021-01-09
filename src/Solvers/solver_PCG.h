

#ifndef UPA_SOLVER_PCG_H
#define UPA_SOLVER_PCG_H


#include "sparse_CSR.h"
#include "solver.h"
#include "preconditioner.h"

namespace upa {

    class Solver_PCG: public Solver {

    private:

        Preconditioner* P;

    public:

        /// Constructors & Destructors
        Solver_PCG(int problemSize, Sparse_CSR* problemMatrix, double *problemVector, Preconditioner *precond);
        ~Solver_PCG() override = default;


        /// Solve the system using the conjugate gradient method
        void solve(double *x0) override;

    };

}


#endif //UPA_SOLVER_PCG_H
