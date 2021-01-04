

#ifndef UPA_SOLVER_AMG_H
#define UPA_SOLVER_AMG_H

#include "sparse_CSR.h"
#include "solver.h"

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
        class Iterator;           // General iterator class
        enum class Iterator_Type; // Enum Class for different iterators
        class Iterator_Jacobi;    // Jacobi iterator
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


    /** General Iterator class
     *
     *  Each iterator has two main functions:
     *      1) A constructor, which sets up the internal variables.
     *      2) An iteration function, which should have no further setup required.
     */
    class Solver_AMG::Iterator {

    public:
        virtual ~Iterator() = default;
        virtual void iterate(double* x_in, double* b, double* x_out) = 0;

    };

    /** Enum Class for different iterators **/
    enum class Solver_AMG::Iterator_Type {
        Jacobi
    };

    /**  Jacobi iterator for AMG
     *   Each iteration performs the action  x <- (1-w) x + w[x + D^-1 (b - Ax)]
     *   where
     *      1) w is a relaxation parameter (tuneable).
     *      2) A is the matrix system, D is its diagonal and b is the system vector.
     */
    class Solver_AMG::Iterator_Jacobi: public Iterator {

    public:
        Iterator_Jacobi(Sparse_CSR* problemMatrix, double relaxationParameter);
        ~Iterator_Jacobi() override = default;
        void iterate(double* x_in, double* b, double* x_out) override;

    private:
        int n;
        Sparse_CSR* A;
        double* D;
        double w;

    };


}

#endif //UPA_SOLVER_AMG_H
