

#ifndef UPA_SOLVER_AMG_H
#define UPA_SOLVER_AMG_H

#include "sparse_CSR.h"
#include "solver.h"
#include "iterator.h"

namespace upa {

    // TODO: Single-level implemented. To have multiple levels, an interpolator has to be implemented (Rugge-Steuben,...)
    //       For now, AMG allows for simple solvers such as Jacobi, Gauss-Seidel, SSOR, ... which are 1-level smoothers.

    /**
     *  Algebraic MultiGrid (AMG) solver for sparse matrices
     **/
    class Solver_AMG: public Solver {

    public:

        /// Constructors & Destructors
        Solver_AMG(int problemSize, Sparse_CSR* problemMatrix, double *problemVector);
        ~Solver_AMG() override = default;

        /// Subclass declaration
        class AMG_Level;          // Level structure for the AMG solver

        /// Additional setup
        void setIteratorType(Iterator_Type iType);

        /// Solve the system using the conjugate gradient method
        void solve(double *x0) override;

    private:
        Iterator_Type _type;

    };

    /** AMG level class
     *
     */
    class Solver_AMG::AMG_Level {
    public:
    private:
    };


}

#endif //UPA_SOLVER_AMG_H
