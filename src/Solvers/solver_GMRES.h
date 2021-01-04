

#ifndef UPA_SOLVER_GMRES_H
#define UPA_SOLVER_GMRES_H

#include "sparse_CSR.h"
#include "solver.h"

namespace upa {

    class Solver_GMRES: public Solver {

    private:

    public:

        /// Constructors & Destructors
        Solver_GMRES(int problemSize, Sparse_CSR* problemMatrix, double *problemVector);
        ~Solver_GMRES() override = default;

        /// Solve the system using the GMRES method
        void solve(double *x0) override;

        void setRestart(bool restart, int m);


    private:

        bool _restart; // Are we using restart?
        int _m;       // Number of steps before restart

        /**@brief Householder matrix calculation
         *
         * Given x[nx], returns w[nx], beta such that
         *
         *    1) w[j0] = 1.0 and w[i] = 0.0 for i < j0
         *    2) (P x)[i] = 0.0 for i > j0 , where P = Id - beta * w·w'
         *
         * @param nx    - Input: Size of x
         * @param j0    - Input: Index where to start orthogonalization
         * @param x     - Input: Input vector
         * @param w     - Output: Houselder vector
         * @param beta  - Output: Householder scalar
         **/
        void _householder(int nx, int j0, double* x, double* w, double& beta);

    };

}

#endif //UPA_SOLVER_CG_H
