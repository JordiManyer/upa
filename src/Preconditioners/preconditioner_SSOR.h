

#ifndef UPA_PRECONDITIONER_SSOR_H
#define UPA_PRECONDITIONER_SSOR_H

#include "preconditioner.h"
#include "sparse_CSR.h"

namespace upa {

    class Preconditioner_SSOR : public Preconditioner {

    public:
        Preconditioner_SSOR(Sparse_CSR *A, double w = 1.0);
        ~Preconditioner_SSOR() override = default;

        void apply(double *r_in, double *r_out) override;

    private:
        int n;
        Sparse_CSR* M;

    };

}

#endif //UPA_PRECONDITIONER_SSOR_H