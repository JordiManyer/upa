
#ifndef UPA_PRECONDITIONER_H
#define UPA_PRECONDITIONER_H

#include <string>

namespace upa {

    class preconditioner {
    private:
        std::string type;
        int n;
        double *D;
        sparse_CSR *M;

    public:
        // Constructor
        preconditioner();

        // Preconditioners
        void jacobi(sparse_CSR *A);           // Jacobi Preconditioner
        void ssor(sparse_CSR *A, double w);  // Symmetric SOR Preconditioner

        // Solve Preconditioner
        void solve(double *x, double *y);         // Generic
        void solve_jacobi(double *x, double *y);  // Jacobi
        void solve_ssor(double *x, double *y);    // SSOR
    };

}

#endif //UPA_PRECONDITIONER_H
