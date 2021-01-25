
#ifndef UPA_REFELEM_NEDELEC_H
#define UPA_REFELEM_NEDELEC_H

#include "elementDefinitions.h"
#include "referenceElement.h"

namespace upa {

    template<>
    class RefElem<ElemType::Triangle, BFType::Nedelec, 1> : public ReferenceElement {
    public:
        RefElem();
        ~RefElem() override = default;

        void evaluateBFs(const double* coords, double *bf) override;
        void evaluateDBFs(const double* coords, double *dbf) override;

    };


}

#endif //UPA_REFELEM_NEDELEC_H
