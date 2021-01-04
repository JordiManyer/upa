
#include "solver_AMG.h"
#include <cmath>
#include <stdexcept>

namespace upa {

/*****************************************************************************************
 ******************              CONSTRUCTORS & DESTRUCTOR             *******************
 *****************************************************************************************/

    Solver_AMG::Solver_AMG(int problemSize, Sparse_CSR* problemMatrix, double *problemVector) {

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
        _type = Iterator_Type::Jacobi;
    }


/*****************************************************************************************
 *************************        OPTIONS AND PARAMETERS         *************************
 *****************************************************************************************/

    void Solver_AMG::setIteratorType(Iterator_Type iType) {
        _type = iType;
    }


/*****************************************************************************************
 *************************                 SOLVER                *************************
 *****************************************************************************************/


    void Solver_AMG::solve(double *x0) {

        /// Initialisation

        /// Main Loop :

        /// Save outputs:

    }


    ///*************************** ITERATORS ***************************///

    Solver_AMG::Iterator_Jacobi::Iterator_Jacobi(Sparse_CSR* problemMatrix, double relaxationParameter) {
        n = problemMatrix->n;
        A = problemMatrix;
        w = relaxationParameter;
        D = new double[n];
        A->getDiag(D);
    }

    void Solver_AMG::Iterator_Jacobi::iterate(double* x_in, double* b, double* x_out) {
        A->matvec(x_in, x_out);
        for (int i = 0; i < n; ++i) {
            x_out[i] = (b[i] - x_out[i])/D[i] + x_in[i];
            x_out[i] = w * x_out[i] + (1.0 - w) * x_in[i];
        }
    }


}