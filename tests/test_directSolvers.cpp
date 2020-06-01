using namespace std;

#include <iostream>
#include "solvers/directSolvers.h"
#include "solvers/conjugateGradient.h"


void testLU();
void testChol();
void testCG();

int main() {
    cout << "Testing Cholesky: " << endl;
    testChol();

    cout << endl << "Testing LU: " << endl;
    testLU();

    cout << endl << "Testing CG: " << endl;
    testCG();

    return 0;
}



void testChol() {
    int n = 3;
    double M[n*n], L[n*n], b[n], x[n];

    // Input
    M[0] = 2.0; M[1] = -1.0; M[2] = 0.0;
    M[3] = -1.0; M[4] = 2.0; M[5] = -1.0;
    M[6] = 0.0; M[7] = -1.0; M[8] = 2.0;

    b[0] = 1.0; b[1] = 2.0; b[2] = 3.0;


    // Solve system
    Chol(n, M, L);
    solveChol(n, L, b, x);


    // Output
    cout << "L = " << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << L[i*n + j] << " , ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "x = " << endl;
    for (int i = 0; i < n; ++i) {
        cout << x[i] << " , ";
    }
    cout << endl;

    double err = 0; double aux;
    for (int i = 0; i < n; ++i) {
        aux = 0;
        for (int j = 0; j < n; ++j) {
            aux += M[i*n+j]*x[j];
        }
        err += (b[i]-aux)*(b[i]-aux);
    }
    cout << "Error : " << err << endl;
}

void testLU() {
    int n = 3;
    double M[n*n], b[n], x[n];
    int P[n];

    // Input
    M[0] = 2.0; M[1] = -1.0; M[2] = 0.0;
    M[3] = -1.0; M[4] = 2.0; M[5] = -1.0;
    M[6] = 0.0; M[7] = -1.0; M[8] = 2.0;

    b[0] = 1.0; b[1] = 2.0; b[2] = 3.0;


    // Solve system
    LU(n,M,P);
    solveLU(n, M , P , b, x);


    // Output
    cout << "L = " << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << M[i*n + j] << " , ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "P = " << endl;
    for (int i = 0; i < n; ++i) {
        cout << P[i] << " , ";
    }
    cout << endl;

    cout << "x = " << endl;
    for (int i = 0; i < n; ++i) {
        cout << x[i] << " , ";
    }
    cout << endl;
}

void testCG() {
    int n = 3;
    double M[n*n], b[n], x[n], x0[n], r[n];

    // Input
    M[0] = 2.0; M[1] = -1.0; M[2] = 0.0;
    M[3] = -1.0; M[4] = 2.0; M[5] = -1.0;
    M[6] = 0.0; M[7] = -1.0; M[8] = 2.0;

    b[0] = 1.0; b[1] = 2.0; b[2] = 3.0;

    x0[0] = 1.0; x0[1] = 1.0; x0[2] = 1.0;


    // Solve system
    conjugateGradient CG(n, M, b);
    CG.configure(0.000001, false, 1);
    CG.solve(x0);

    CG.getSolution(x);
    CG.getResidual(r);

    double err = 0; double aux;
    for (int i = 0; i < n; ++i) {
        aux = 0;
        for (int j = 0; j < n; ++j) {
            aux += M[i*n+j]*x[j];
        }
        err += (b[i]-aux)*(b[i]-aux);
    }
    cout << "Error1 : " << err << endl;

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
