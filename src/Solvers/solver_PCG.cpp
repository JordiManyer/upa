

#include "solver_PCG.h"
#include <cmath>


namespace upa {

/*****************************************************************************************
 ******************              CONSTRUCTORS & DESTRUCTOR             *******************
 *****************************************************************************************/

    Solver_PCG::Solver_PCG(int problemSize, Sparse_CSR* problemMatrix, double *problemVector, Preconditioner *precond) {

        // Inputs
        n = problemSize;
        A = problemMatrix;
        b = problemVector;
        P = precond;

        // Outputs
        hasConverged = false;
        sol = new double[n];
        residual = new double[n];
        numIter = -1;
        finalError = -1;

        // Set default parameters
        tol = 1.e-12;
        maxIter = 2 * n;
        verbose = 0;
    }


/*****************************************************************************************
 *************************        OPTIONS AND PARAMETERS         *************************
 *****************************************************************************************/



/*****************************************************************************************
 *************************       CONJUGATE GRADIENT METHOD       *************************
 *****************************************************************************************/


    void Solver_PCG::solve(double *x0) {
        int k;
        double rk[n], rkold[n], dk[n], Adk[n], xk[n], zk[n];
        double alpha, beta, err, aux;

        /// Initialisation
        k = 0;
        err = 0.0;
        A->matvec(x0, Adk);
        for (int i = 0; i < n; ++i) {
            xk[i] = x0[i];
            rk[i] = b[i] - Adk[i];
            err += rk[i] * rk[i];
        }
        P->apply(rk,zk); // z0 = P·r0
        for (int i = 0; i < n; ++i) dk[i] = zk[i]; // d0 = z0
        err = sqrt(err);

        /// Main Loop :
        while (err > tol and k < maxIter) {
            // Calculate A * dk
            A->matvec(dk, Adk);

            // alpha = (zk·rk) / (dk·A·dk)
            alpha = 0.0;
            aux = 0;
            for (int i = 0; i < n; ++i) {
                alpha += zk[i] * rk[i];
                aux += dk[i] * Adk[i];
            }
            alpha /= aux;

            // xkp1 = xk + alpha * dk ; rkp1 = rk - alpha * A * dk
            for (int i = 0; i < n; ++i) {
                xk[i] += alpha * dk[i];

                rkold[i] = rk[i];
                rk[i] -= alpha * Adk[i];
                Adk[i] = zk[i]; // saving old zk
            }

            // Update stop conditions
            err = 0.0;
            for (int i = 0; i < n; ++i) {
                err += rk[i] * rk[i];
            }
            if (err < tol) break;

            // zkp1 = P·zk
            P->apply(rk,zk);

            // beta = (zkp1·rkp1) / (zk·rk)
            beta = 0.0;
            aux = 0.0;
            for (int i = 0; i < n; ++i) {
                beta += zk[i] * rk[i];
                aux += Adk[i] * rkold[i];
            }
            beta /= aux;

            // dkp1 = zkp1 + beta * dk
            for (int i = 0; i < n; ++i) dk[i] = zk[i] + beta * dk[i];

            ++k;
        }

        /// Save outputs:
        if (k == n * n) hasConverged = false;
        else hasConverged = true;

        for (int i = 0; i < n; ++i) {
            sol[i] = xk[i];
            residual[i] = rk[i];
        }
        numIter = k;
        finalError = err;
    }


}