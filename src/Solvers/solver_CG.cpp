
#include "solver_CG.h"
#include <cmath>

namespace upa {

/*****************************************************************************************
 ******************              CONSTRUCTORS & DESTRUCTOR             *******************
 *****************************************************************************************/

    Solver_CG::Solver_CG(int problemSize, Sparse_CSR* problemMatrix, double *problemVector) {

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
        maxIter = 2 * n;
        verbose = 0;
    }


/*****************************************************************************************
 *************************        OPTIONS AND PARAMETERS         *************************
 *****************************************************************************************/



/*****************************************************************************************
 *************************       CONJUGATE GRADIENT METHOD       *************************
 *****************************************************************************************/


    void Solver_CG::solve(double *x0) {
        int k;
        double rk[n], rkold[n], dk[n], Adk[n], xk[n];
        double alpha, beta, err, aux;

        /// Initialisation
        k = 0;
        err = 0.0;
        A->matvec(x0, Adk);
        for (int i = 0; i < n; ++i) {
            xk[i] = x0[i];
            rk[i] = b[i] - Adk[i];
            dk[i] = rk[i];
            err += rk[i] * rk[i];
        }
        err = sqrt(err);

        /// Main Loop :
        while (err > tol and k < maxIter) {
            // Calculate A * dk
            A->matvec(dk, Adk);

            // alpha = (dk·rk) / (dk·A·dk)
            alpha = 0.0;
            aux = 0;
            for (int i = 0; i < n; ++i) {
                alpha += dk[i] * rk[i];
                aux += dk[i] * Adk[i];
            }
            alpha /= aux;

            // xkp1 = xk + alpha * dk ; rkp1 = rk - alpha * A * dk
            for (int i = 0; i < n; ++i) {
                xk[i] += alpha * dk[i];

                rkold[i] = rk[i];
                rk[i] -= alpha * Adk[i];
            }

            // beta = (rkp1·rkp1) / (rk·rk)
            beta = 0.0;
            aux = 0.0;
            for (int i = 0; i < n; ++i) {
                beta += rk[i] * rk[i];
                aux += rkold[i] * rkold[i];
            }
            beta /= aux;

            // dkp1 = rkp1 + beta * dk
            for (int i = 0; i < n; ++i) dk[i] = rk[i] + beta * dk[i];

            // Update stop conditions
            err = 0.0;
            double orthogonal = 0;
            for (int i = 0; i < n; ++i) {
                orthogonal += rk[i] * rkold[i];
                err += rk[i] * rk[i];
            }
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