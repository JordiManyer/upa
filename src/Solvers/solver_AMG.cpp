
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

}