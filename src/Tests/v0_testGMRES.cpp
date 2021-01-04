
#include <iostream>
#include <cmath>
#include "sparse_LIL.h"
#include "sparse_CSR.h"
#include "solver_GMRES.h"
#include "debugIO.h"

using namespace std;
using namespace upa;


int main() {
    int n = 10;
    Sparse_LIL lil(n);

    ///  Assembling the matrix in LIL format
    for (int i = 0; i < n; ++i) {
        lil.put(i,i,2);
        if (i != 0  ) lil.put(i,i-1,-1);
        if (i != n-1) lil.put(i,i+1,-1);
    }

    lil.put(0,2,3); // Break symmetry

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
        b[i] = i;
        x0[i] = 1.0;
    }

    // Solve system
    Solver_GMRES GMRES(n, &csr, b);
    GMRES.setTolerance(0.000001);
    GMRES.setRestart(true, 5);
    GMRES.setVerbosity(1);
    GMRES.solve(x0);

    GMRES.getSolution(x);
    GMRES.getResidual(r);

    // Output
    cout << "Converged : " << GMRES.getConvergence() << endl;
    cout << "Number of iterations: " << GMRES.getNumIter() << endl;
    cout << "Error2 : " << GMRES.getError() << endl;

    cout << "x = " << endl;
    printArray(n,x);

    // Check that the residual is well calculated
    double err = 0.0;
    csr.matvec(x,r);
    for (int i = 0; i < n; ++i) {
        r[i] -= b[i];
        err += r[i]*r[i];
    }
    err = sqrt(err);
    cout << "Residual = " << err << endl;
    printArray(n,r);
}