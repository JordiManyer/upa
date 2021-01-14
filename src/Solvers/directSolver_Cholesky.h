

#ifndef UPA_DIRECTSOLVER_CHOLESKY_H
#define UPA_DIRECTSOLVER_CHOLESKY_H

#include "directSolver.h"
#include <stack>

namespace upa {

    class DirectSolver_Cholesky : public DirectSolver {

    public:
        DirectSolver_Cholesky(int problemSize, Sparse_CSR* problemMatrix, double *problemVector);
        ~DirectSolver_Cholesky() override = default;

        void solve(double *x0) override;

        // Elimination tree
        int *parent;
        std::vector<int>* descendants;
        int *postorder;
        int *postorder_inv;

        // Skeleton Graph
        std::vector<int> *leaves;

        // Cholesky matrix L
        int nnz;
        int *rows;
        int *cols;
        double *values;


        /** @brief Produces the elimination tree of the CSR matrix A.
         *
         *  The tree is saved in the array 'parent', which at the end
         *  of the algorithm holds at each position parent[i] the index
         *  of its parent node in the elimination tree.
         *
         *  Follows Algorithm 4.1 in [Liu,1986]
         */
        void _getEliminationTree();

        /** @brief Generates a postorder ordering of the elimination tree.
         *
         *  Uses a DFS-like algorithm to avoid recursion.
         *  See [Aho, Data Structures and Algorithms] and [Duff,1983].
         */
        void _postOrder();

        /** @brief Selects leaf nodes for the skeleton graph.
         *
         *  Uses Corollary 4.3 in [Liu,1986].
         */
        void _findLeaves();

        /** @brief Creates sparsity pattern of the Cholesky matrix L
         *
         */
        void _symbolicFactorization();

        void _factorize();

    };

}

#endif //UPA_DIRECTSOLVER_CHOLESKY_H
