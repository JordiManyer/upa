

#ifndef UPA_PRECONDITIONER_DIAG_H
#define UPA_PRECONDITIONER_DIAG_H

#include "preconditioner.h"
#include "sparse_CSR.h"

namespace upa {

    class Preconditioner_Diag : public Preconditioner {

    public:
        Preconditioner_Diag(Sparse_CSR *A);

        ~Preconditioner_Diag() override = default;

        void apply(double *r_in, double *r_out) override;

    private:
        int n;
        double *D;

    };

}

#endif //UPA_PRECONDITIONER_DIAG_H
