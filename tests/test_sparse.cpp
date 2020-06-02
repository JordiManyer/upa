using namespace std;

#include <iostream>
#include "solvers/directSolvers.h"
#include "solvers/conjugateGradient.h"
#include "sparse/sparse.h"

void testSparse();
void testSparseCG();

int main() {

    testSparse();
    testSparseCG();

    return 0;
}


void testSparse() {
    int n = 4;
    sparse_LIL lil(n);

    ///  Assembling the matrix
    lil.assemble_elem(1, 2, 0.5);
    lil.assemble_elem(2, 2, -0.5);
    lil.assemble_elem(3, 1, 1.0);
    lil.assemble_elem(1, 2, 0.5); // The repeated element will be summed to the other contribution
    lil.assemble_elem(0, 3, -2.0);
    lil.assemble_elem(0, 0, 2.5);

    /// Printing LIL matrix in 2 ways
    cout << "List of elements: " << endl;
    for (int i = 0; i < n; ++i) {
        cout << i << "->   ";
        for (auto it = lil.A[i].begin() ; it < lil.A[i].end() ; ++it) {
            cout << " (" << it->first << ';' << it->second << ") ," ;
        }
        cout << endl;
    }
    cout << endl;

    cout << "LIL Complete matrix: " << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n ; ++j) {
            cout << "    " << lil.get_elem_by_index(i,j);
        }
        cout << endl;
    }
    cout << endl;

    /// From LIL to CSR
    sparse_CSR* csr =  LILtoCSR(&lil);

    /// Printing CSR matrix in 2 ways
    cout << "CSR compact form: " << endl;
    for (int i = 0; i < csr->n+1; ++i) {
        cout << csr->rows[i] << ' ' ;
    }
    cout << endl;
    for (int i = 0; i < csr->nnz; ++i) {
        cout << csr->cols[i] << ' ' ;
    }
    cout << endl;
    for (int i = 0; i < csr->nnz; ++i) {
        cout << csr->values[i] << ' ' ;
    }
    cout << endl;

    cout << "CSR Complete matrix: " << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n ; ++j) {
            cout << "    " << csr->get(i,j);
        }
        cout << endl;
    }
    cout << endl;
}


void testSparseCG() {
    int n = 3;
    sparse_LIL lil(n);

    ///  Assembling the matrix
    for (int i = 0; i < n; ++i) {
        lil.assemble_elem(i,i,2);
        if (i != 0  ) lil.assemble_elem(i,i-1,-1);
        if (i != n-1) lil.assemble_elem(i,i+1,-1);
    }

    /// From LIL to CSR
    sparse_CSR* csr =  LILtoCSR(&lil);
    cout << "Complete matrix: " << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n ; ++j) {
            cout << "    " << csr->get(i,j);
        }
        cout << endl;
    }
    cout << endl;

    // System declaration
    double b[n], x[n], x0[n], r[n];
    b[0] = 1.0; b[1] = 2.0; b[2] = 3.0;
    x0[0] = 1.0; x0[1] = 1.0; x0[2] = 1.0;

    // Solve system
    conjugateGradient CG(n, csr, b);
    CG.configure(0.000001, true, 1);
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



