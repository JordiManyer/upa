
#include "solver_GMRES.h"
#include <cmath>

namespace upa {

/*****************************************************************************************
 ******************              CONSTRUCTORS & DESTRUCTOR             *******************
 *****************************************************************************************/

    Solver_GMRES::Solver_GMRES(int problemSize, Sparse_CSR* problemMatrix, double *problemVector) {

        // Inputs
        n = problemSize;
        A = problemMatrix;
        b = problemVector;

        // Outputs
        hasConverged = false;
        sol = new double[n];
        residual = new double[n];
        numIter = -1;
        finalError = -1;

        // Set default parameters
        tol = 1.e-12;
        verbose = 0;
        _m = n;
        _restart = false;
    }


/*****************************************************************************************
 *************************        OPTIONS AND PARAMETERS         *************************
 *****************************************************************************************/

    void Solver_GMRES::setRestart(bool restart, int m) {
        _m = m;
        _restart = restart;
        if (m > n) _m = n;
        if (not restart) _m = n;
    }


/*****************************************************************************************
 *************************          GMRES SOLVER METHOD          *************************
 *****************************************************************************************/


    void Solver_GMRES::solve(double *x0) {
        int m;
        double vj[n], wj[n], xk[n], vAux[n];
        double H[n*(_m+1)], betas[_m+1];
        double beta, rNorm, s;

        /// Initialisation
        m = _m;
        rNorm = 0.0;
        A->matvec(x0, vAux);
        for (int i = 0; i < n; ++i) {
            xk[i] = x0[i];         // Copy initial condition
            vj[i] = b[i] - vAux[i]; // v0 = r0 (initial residual)
            rNorm += vj[i] * vj[i]; // norm(r0)^2
        }
        rNorm = sqrt(rNorm);

        /// Householder loop: Building the orthonormal basis H (QR-factorized)
        for (int j = 0; j < m+1; ++j) {

            // Householder vector of vj
            _householder(n,j,vj,wj,beta);

            // Intermediate variable s = beta * wj' * vj
            s = 1.0 * vj[j];
            for (int i = j+1; i < n; ++i) s += wj[i] * vj[i];
            s *= beta;

            // Calculate and store the jth column of H
            betas[j] = beta;
            for (int i = 0; i < n; ++i) {
                if (i > j) H[i+n*j] = wj[i];
                else H[i+n*j] = vj[i] - s * wj[i]; // hj = Pj * vj = vj - beta * wj * wj' * vj
            }

            // Compute new vj
            for (int i = 0; i < n; ++i) vAux[i] = 0.0; vAux[j] = 1.0; // Init: z = ej
            for (int k = j; k >= 0; --k) { // z = P1 * P2 * ... * Pj * ej
                s = 1.0 * vAux[k];
                for (int i = k+1; i < n; ++i) s += H[i+n*k] * vAux[i];
                s *= betas[k];
                vAux[k] -= s * 1.0;
                for (int i = k+1; i < n; ++i) vAux[i] -= s * H[i+n*k];
            }

            A->matvec(vAux, vj); // z' = A z
            for (int k = 0; k < j+1; ++k) { // vj = Pj * ... * P1 * z'
                s = 1.0 * vj[k];
                for (int i = k+1; i < n; ++i) s += H[i+n*k] * vj[i];
                s *= betas[k];
                vj[k] -= s * 1.0;
                for (int i = k+1; i < n; ++i) vj[i] -= s * H[i+n*k];
            }

        }

        // TODO : Solve the least squares problem using the QR factorization
        // TODO : Restart scheme

    }

    /// ---------------------------------------------------------------------------------------------------------

    void Solver_GMRES::_householder(int nx, int j0, double* x, double* wj, double& beta) {
        double s; double mu;

        s = 0.0;
        wj[j0] = 1.0;
        for (int i = 0; i < j0; ++i) wj[i] = 0.0;
        for (int i = j0+1; i < nx; ++i) {
            s += x[i]*x[i];
            wj[i] = x[i];
        }

        if (s < 1.e-14) beta = 0.0;
        else {
            mu = sqrt(x[j0]*x[j0] + s);
            if (x[j0] <= 0) wj[j0] = x[j0] - mu;
            else wj[j0] = -s / (x[j0] + mu);

            beta = 2.0*wj[j0]*wj[j0]/(s + wj[j0]*wj[j0]);
            for (int i = j0+1; i < nx; ++i) {
                wj[i] /= wj[j0];
            }
            wj[j0] = 1.0;
        }
    }

}