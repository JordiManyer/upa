
#include "directSolver_Cholesky.h"

#include <iostream>
#include "debugIO.h"
using namespace std;

namespace upa {

    DirectSolver_Cholesky::DirectSolver_Cholesky(int problemSize, Sparse_CSR* problemMatrix, double *problemVector) {
        n = problemSize;
        A = problemMatrix;
        b = problemVector;

        // Create elimination tree
        parent = new int[n];
        _getEliminationTree();

        // Find descendants
        descendants = new std::vector<int>[n];
        for (int i = n-2; i >= 0; --i) descendants[parent[i]].push_back(i);

        // Find a post2pre ordering for the elimination tree
        post2pre = new int[n];
        pre2post = new int[n];
        _postOrder();

        // Create skeleton graph
        leaves = new std::vector<int> [n];
        _findLeaves();

        // Creates the structure for L
        _findSkeleton();

        // Fill the cholesky matrix L
        _factorize();

    }

    void DirectSolver_Cholesky::solve(double *x0) {

    }


    /// TODO: Compression can enhance algorithm, see [TARJAN, Data Structures and Network Algorithms. 1983.]
    void DirectSolver_Cholesky::_getEliminationTree() {
        int r, temp;
        int ancestor[n];

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


    void DirectSolver_Cholesky::_postOrder() {
        int k; int pos = n-1;
        stack<int> s; s.push(n-1);
        while (not s.empty()) {
            k = s.top();
            s.pop();
            for (int & i : descendants[k]) s.push(i);
            post2pre[pos] = k;
            pre2post[k] = pos;
            --pos;
        }
    }


    void DirectSolver_Cholesky::_findLeaves() {
        // Find neighbors and order in post2pre
        for (int i = 0; i < n; ++i) {
            std::vector<int> neighbors;
            for (int j = A->rows[post2pre[i]]; j < A->rows[post2pre[i]+1]; ++j) {
                if(pre2post[A->cols[j]] < i) neighbors.push_back(pre2post[A->cols[j]]);
            }
            sort(neighbors.begin(), neighbors.end());

            // Select direct ancestors using Corollary 4.3 in [Liu,1986].
            for (int j = 0; j < neighbors.size(); ++j) {
                if (j == 0 or neighbors[j-1] < neighbors[j] - descendants[post2pre[neighbors[j]]].size()) leaves[i].push_back(neighbors[j]);
            }
        }
    }


    void DirectSolver_Cholesky::_findSkeleton() {
        int d, dNext, k;

        // First traverse of the elimination tree: Count nonzero entries of L
        nnz = 0;
        rows = new int [n+1]; rows[0] = 0;
        for (int i = 0; i < n; ++i) {
            rows[i+1] = 0;

            for (auto ptr = leaves[i].begin(); ptr != leaves[i].end(); ++ptr) {
                d = *ptr;

                if (std::next(ptr,1) != leaves[i].end()) dNext = *std::next(ptr,1);
                else dNext = i;

                while (d < dNext) {
                    rows[i+1]++;
                    nnz++;
                    d = pre2post[parent[post2pre[d]]];
                }
            }
        }
        for (int i = 0; i < n; ++i) rows[i+1] += rows[i]; // Cumulated counts

        // Allocate L
        cols = new int[nnz];
        values = new double[nnz];

        // Second traverse of the elimination tree: find nonzero indexes of L
        for (int i = 0; i < n; ++i) {
            k = rows[i];

            for (auto ptr = leaves[i].begin(); ptr != leaves[i].end(); ++ptr) {
                d = *ptr;

                if (std::next(ptr,1) != leaves[i].end()) dNext = *std::next(ptr,1);
                else dNext = i;

                while (d < dNext) {
                    cols[k] = d;
                    k++;
                    d = pre2post[parent[post2pre[d]]];
                }
            }
        }
    }


    void DirectSolver_Cholesky::_factorize() {

    }


}