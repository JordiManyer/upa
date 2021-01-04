
#include <iostream>
#include "sparse_LIL.h"
#include "sparse_CSR.h"
#include "solver_AMG.h"

using namespace std;
using namespace upa;


int main() {
    int n = 3;
    Sparse_LIL lil(n);

    ///  Assembling the matrix in LIL format
    for (int i = 0; i < n; ++i) {
        lil.put(i,i,2);
        if (i != 0  ) lil.put(i,i-1,-1);
        if (i != n-1) lil.put(i,i+1,-1);
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
    double b[n], xkp1[n], xk[n], r[n];
    b[0] = 1.0; b[1] = 2.0; b[2] = 3.0;
    xk[0] = 1.0; xk[1] = 1.0; xk[2] = 1.0;

    // Solve system
    double tol = 1.e-5;
    int maxIter = 100;

    auto iterator = new Solver_AMG::Iterator_Jacobi(&csr,0.3);

    double error = 10.0; int k = 0;
    while (error > tol and k < maxIter) {

        iterator->iterate(xk,b,xkp1);

        for (int i = 0; i < n; ++i) {
            cout << xkp1[i] << " , ";
        }
        cout << endl;

        error = 0.0;
        for (int i = 0; i < n; ++i) {
            error += (xk[i]-xkp1[i]) * (xk[i]-xkp1[i]);
            xk[i] = xkp1[i];
        }
        k++;
    }

    if (k < maxIter) cout << "The algorithm has converged in " << k << " iterations!" << endl;
    else cout << "The algorithm has NOT converged!" << endl;

    cout << "Final error = " << error << endl;

    cout << "x = " << endl;
    for (int i = 0; i < n; ++i) {
        cout << xkp1[i] << " , ";
    }
    cout << endl;

}