
#include "solver_GMRES.h"
#include <cmath>
#include "debugIO.h"

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
        maxIter = n*n;
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

    /** GMRES Algorithm
     * Implementation follows Algorithm (6.10) from https://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf
     * except for the definition of the Householder vectors. A slight modification has been used that allows
     * compact storing of the information. This modification of the Householder transformations can be found
     * in http://eprints.ma.man.ac.uk/1192/1/covered/MIMS_ep2008_111.pdf
     *
     * @param x0 - Initial approximation
     */
    void Solver_GMRES::solve(double *x0) {
        int m;
        double vj[n], wj[n], xk[n], vAux[n];
        double H[n*(_m+1)], betas[_m+1];
        double g[_m+1];
        double beta, rNorm, s, c;

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
        hasConverged = (rNorm < tol);

        /// Restart Loop: Updates xk and rNorm at each step
        int k = 0;
        while (not hasConverged and k < maxIter) {

            /// Householder loop: Builds the orthonormal basis for the Krylov space
            for (int j = 0; j < m + 1; ++j) {

                // Householder vector of vj
                _householder(n, j, vj, wj, beta);

                // Intermediate variable s = beta * wj' * vj
                s = 1.0 * vj[j];
                for (int i = j + 1; i < n; ++i) s += wj[i] * vj[i];
                s *= beta;

                // Calculate and store the jth column of H
                betas[j] = beta;
                for (int i = 0; i < n; ++i) {
                    if (i > j) H[i + n * j] = wj[i];
                    else H[i + n * j] = vj[i] - s * wj[i]; // hj = Pj * vj = vj - beta * wj * wj' * vj
                }

                // Compute new vj
                for (int i = 0; i < n; ++i) vAux[i] = 0.0;
                vAux[j] = 1.0; // Init: z = ej
                for (int l = j; l >= 0; --l) { // z = P0 * P1 * ... * Pj * ej
                    s = 1.0 * vAux[l];
                    for (int i = l + 1; i < n; ++i) s += H[i + n * l] * vAux[i];
                    s *= betas[l];
                    vAux[l] -= s * 1.0;
                    for (int i = l + 1; i < n; ++i) vAux[i] -= s * H[i + n * l];
                }

                A->matvec(vAux, vj); // z' = A z
                for (int l = 0; l < j + 1; ++l) { // vj = Pj * ... * P0 * z'
                    s = 1.0 * vj[l];
                    for (int i = l + 1; i < n; ++i) s += H[i + n * l] * vj[i];
                    s *= betas[l];
                    vj[l] -= s * 1.0;
                    for (int i = l + 1; i < n; ++i) vj[i] -= s * H[i + n * l];
                }

            } // End of Householder loop

            /// Solve the Least Squares optimization problem y_m = ArgMin Norm[ e_1' h_0 - H_m y ]
            // Init g_m = e_1' h_0 = [h_0]_0
            for (int i = 0; i < m + 1; ++i) g[i] = 0.0;
            g[0] = H[0 + n * 0];

            // Triangularize H_m into R_m by using plane rotations
            for (int i = 0; i < m; ++i) {
                // Setup rotation to eliminate element [H_m]_ii
                beta = sqrt(H[i + n * (i + 1)] * H[i + n * (i + 1)] + H[(i + 1) + n * (i + 1)] * H[(i + 1) + n * (i + 1)]);
                s = H[(i + 1) + n * (i + 1)] / beta;
                c = H[i + n * (i + 1)] / beta;

                // Apply rotation to g_m
                beta = c * g[i] + s * g[i + 1]; // Placeholder for g[i]
                g[i + 1] = -s * g[i] + c * g[i + 1];
                g[i] = beta;

                // Apply rotation to H_m
                for (int j = i + 1; j < m + 1; ++j) {
                    beta = c * H[i + n * j] + s * H[(i + 1) + n * j]; // Placeholder for H[i+n*j]
                    H[(i + 1) + n * j] = -s * H[i + n * j] + c * H[(i + 1) + n * j];
                    H[i + n * j] = beta;
                }
            }

            // Solving triangular system R_m y = g_m
            for (int i = m - 1; i >= 0; --i) {
                for (int j = i + 1; j < m; ++j) g[i] -= H[i + n * (j + 1)] * g[j];
                g[i] /= H[i + n * (i + 1)];
            }

            // Undo the basis change to find the new iterate
            for (int i = 0; i < n; ++i) vAux[i] = 0.0; // Init: z = 0
            for (int j = m - 1; j >= 0; --j) { // z = Pj (y[j] ej + z) for j = m, ..., 1
                vAux[j] += g[j];

                s = 1.0 * vAux[j];
                for (int i = j + 1; i < n; ++i) s += H[i + n * j] * vAux[i];
                s *= betas[j];

                vAux[j] -= s * 1.0;
                for (int i = j + 1; i < n; ++i) vAux[i] -= s * H[i + n * j];
            }
            for (int i = 0; i < n; ++i) xk[i] += vAux[i];

            // Update residual
            rNorm = fabs(g[m]);
            A->matvec(xk, vAux);
            for (int i = 0; i < n; ++i) vj[i] = b[i] - vAux[i];

            // Update stop conditions
            hasConverged = (rNorm < tol);
            k++;
        } // End of restart loop

        /// Copy results
        for (int i = 0; i < n; ++i) sol[i] = xk[i];
        finalError = rNorm;
        numIter = k;
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