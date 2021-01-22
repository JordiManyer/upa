
#ifndef UPA_REFELEM_LAGRANGIAN_H
#define UPA_REFELEM_LAGRANGIAN_H

#include "elementDefinitions.h"
#include "referenceElement.h"

namespace upa {

    template<>
    class RefElem<ElemType::Line, BFType::Lagrangian, 1> : public ReferenceElement {
    public:
        RefElem();
        ~RefElem() override = default;

        void evaluateBFs(const double* coords, double *bf) override;
        void evaluateDBFs(const double* coords, double *dbf) override;

    };

    template<>
    class RefElem<ElemType::Line, BFType::Lagrangian, 2> : public ReferenceElement {
    public:
        RefElem();
        ~RefElem() override = default;

        void evaluateBFs(const double* coords, double *bf) override;
        void evaluateDBFs(const double* coords, double *dbf) override;

    };


    template<>
    class RefElem<ElemType::Square, BFType::Lagrangian, 1> : public ReferenceElement {
    public:
        RefElem();
        ~RefElem() override = default;

        void evaluateBFs(const double* coords, double *bf) override;
        void evaluateDBFs(const double* coords, double *dbf) override;

    };


    template<>
    class RefElem<ElemType::Square, BFType::Lagrangian, 2> : public ReferenceElement {
    public:
        RefElem();
        ~RefElem() override = default;

        void evaluateBFs(const double* coords, double *bf) override;
        void evaluateDBFs(const double* coords, double *dbf) override;

    };



}



#endif //UPA_REFELEM_LAGRANGIAN_H
