

#include "directSolver.h"


namespace upa {

    void DirectSolver::setVerbosity(int verbosity) {
        verbose = verbosity;
    }

    void DirectSolver::getSolution(double *x) {
        for (int i = 0; i < n; ++i) x[i] = sol[i];
    }

}