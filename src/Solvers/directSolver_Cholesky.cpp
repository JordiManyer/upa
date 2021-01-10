
#include "directSolver_Cholesky.h"

#include <iostream>
#include "debugIO.h"
using namespace std;

namespace upa {

    DirectSolver_Cholesky::DirectSolver_Cholesky(int problemSize, Sparse_CSR* problemMatrix, double *problemVector) {
        n = problemSize;
        A = problemMatrix;
        b = problemVector;

        parent = new int[n];
        _getEliminationTree();
    }

    void DirectSolver_Cholesky::solve(double *x0) {

    }


    void DirectSolver_Cholesky::_getEliminationTree() {
        int r, temp;
        int ancestor[n];
        printArray(n+1,A->rows);
        printArray(A->nnz,A->cols);
        for (int i = 0; i < n; ++i) {
            parent[i] = -1;
            ancestor[i] = -1;

            // For x_k in Adj(x_i) and k < i
            for (int k = A->rows[i]; k < A->rows[i+1] and A->cols[k] < i; ++k) {
                // Find the root x_r of the tree in the forest containing the node x_{A->cols[k]}
                r = A->cols[k];
                while (ancestor[r] != -1 and ancestor[r] != i) {
                    temp = ancestor[r];
                    ancestor[r] = i;
                    r = temp;
                }
                if (ancestor[r] == -1) {
                    ancestor[r] = i;
                    parent[r] = i;
                }
            }
        }
    }

    


}