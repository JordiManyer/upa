
#include <iostream>
#include "sparse_LIL.h"
#include "sparse_CSR.h"
#include "solver_CG.h"

using namespace std;
using namespace upa;


int main() {
    int n = 4;
    Sparse_LIL lil(n);

    ///  Assembling the matrix in LIL format
    for (int i = 0; i < n; ++i) {
        lil.put(i,i,10.0);
        if (i != 0  ) lil.put(i,i-1,-5.0);
        if (i != n-1) lil.put(i,i+1,-5.0);
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
        b[i] = 0.2;
        x0[i] = 1.0;
    }

    // Solve system
    Solver_CG CG(n, &csr, b);
    CG.setTolerance(0.000001);
    CG.setVerbosity(1);
    CG.solve(x0);

    CG.getSolution(x);
    CG.getResidual(r);

    // Output
    cout << "Converged : " << CG.getConvergence() << endl;
    cout << "Number of iterations: " << CG.getNumIter() << endl;
    cout << "Error2 : " << CG.getError() << endl;

    cout << "x = " << endl;
    for (int i = 0; i < n; ++i) {
        cout << x[i] << " , ";
    }
    cout << endl;

}