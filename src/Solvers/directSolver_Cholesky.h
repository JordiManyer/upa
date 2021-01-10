

#ifndef UPA_DIRECTSOLVER_CHOLESKY_H
#define UPA_DIRECTSOLVER_CHOLESKY_H

#include "directSolver.h"

namespace upa {

    class DirectSolver_Cholesky : public DirectSolver {

    public:
        DirectSolver_Cholesky(int problemSize, Sparse_CSR* problemMatrix, double *problemVector);
        ~DirectSolver_Cholesky() override = default;

        void solve(double *x0) override;


        int *parent;

        /** @brief Produces the elimination tree of the CSR matrix A.
         **
         **  The tree is saved in the array 'parent', which at the end
         **  of the algorithm holds at each position parent[i] the index
         **  of its parent node in the elimination tree.
         **
         **  Follows Algorithm 4.1 in [Liu,1986]
         */
        void _getEliminationTree();


    };

}

#endif //UPA_DIRECTSOLVER_CHOLESKY_H
