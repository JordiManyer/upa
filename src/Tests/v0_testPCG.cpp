

#include <iostream>
#include "sparse_LIL.h"
#include "sparse_CSR.h"
#include "solver_PCG.h"
#include "preconditioner_diag.h"
#include "debugIO.h"

using namespace std;
using namespace upa;


int main() {
    int n = 5;
    Sparse_LIL lil(n);

    ///  Assembling the matrix in LIL format
    for (int i = 0; i < n; ++i) {
        lil.put(i,i,2.0/n);
        if (i != 0  ) lil.put(i,i-1,-1.0/n);
        if (i != n-1) lil.put(i,i+1,-1.0/n);
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

    // Solve system
    Preconditioner_Diag P(&csr);

    Solver_PCG PCG(n, &csr, b, &P);
    PCG.setTolerance(0.000001);
    PCG.setVerbosity(1);
    PCG.solve(x0);

    PCG.getSolution(x);
    PCG.getResidual(r);

    // Output
    cout << "Converged : " << PCG.getConvergence() << endl;
    cout << "Number of iterations: " << PCG.getNumIter() << endl;
    cout << "Error2 : " << PCG.getError() << endl;

    cout << "x = " << endl;
    printArray(n,x);

}