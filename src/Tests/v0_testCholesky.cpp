
#include <iostream>
#include "sparse_LIL.h"
#include "sparse_CSR.h"
#include "directSolver_Cholesky.h"
#include "debugIO.h"
#include "denseSolvers.h"

using namespace std;
using namespace upa;


int main() {
    int n = 6;
    Sparse_LIL lil(n);

    ///  Assembling the matrix in LIL format
/*    for (int i = 0; i < n; ++i) {
        lil.put(i,i,2.0);
        if (i != 0  ) lil.put(i,i-1,-1.0);
        if (i != n-1) lil.put(i,i+1,-1.0);
    }*/

    double dense[36] = {2, 0, 1, 0, 1, 0,
                        0, 2, 0, 0, 1, 1,
                        1, 0, 2, 1, 1, 0,
                        0, 0, 1, 2, 0, 1,
                        1, 1, 1, 0, 2, 0,
                        0, 1, 0, 1, 0, 2};
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (dense[i*n + j] > 0.1) lil.put(i,j,dense[i*n+j]);
        }
    }

    /// Convert to CSR
    auto csr = Sparse_CSR(&lil);
    cout << "Complete matrix: " << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n ; ++j) {
            cout << "    " << csr(i,j);
        }
        cout << endl;
    }
    cout << endl;

    // System declaration
    double b[n], x[n], x0[n], r[n];
    for (int i = 0; i < n; ++i) {
        b[i] = 1.0;
        x0[i] = 1.0;
    }

    auto chol = new DirectSolver_Cholesky(n,&csr,b);

    cout << "parent = " << endl;
    printArray(n,chol->parent);
    cout << "postorder = " << endl;
    printArray(n,chol->postorder);
    cout << "postorder_inv = " << endl;
    printArray(n,chol->postorder_inv);
    cout << "leaves = " << endl;
    for (int i = 0; i < n; ++i) {
        if (chol->leaves[i].size() == 0) cout << " Empty ";
        for (auto j: chol->leaves[i]) cout << j << "  ";
        cout << endl;
    }
    cout << endl;

    cout << "nnz counts: " << endl;
    printArray(n+1,chol->rows);

    cout << "L skeleton: " << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = chol->rows[i]; j < chol->rows[i+1]; ++j) cout << "   " << chol->cols[j];
    cout << endl;
    }
    cout << endl;


    // Check with dense solver:
    double denseMat[n*n], denseL[n*n];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            denseMat[i*n+j] = lil(i,j);
        }
    }
    Chol(n,denseMat,denseL);
    printMatrix(n,n,denseL);


}