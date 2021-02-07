
#ifndef UPA_REFELEM_NEDELEC_H
#define UPA_REFELEM_NEDELEC_H

#include "elementDefinitions.h"
#include "referenceElement.h"
#include "refElem_Lagrangian.h"

namespace upa {

    template<>
    class RefElem<ElemType::Triangle, BFType::Nedelec, 1> : public ReferenceElement {
    public:
        RefElem();
        ~RefElem() override = default;

        void evaluateBFs(const double* refCoords, double *bf) override;
        void evaluateDBFs(const double* refCoords, double *dbf) override;

        void getJacobian(int iG, const double* nodeCoords, double* J) override;
        void getJacobian(const double* dbf, const double* nodeCoords, double* J) override;

        void getPhysicalCoords(int iG, const double* nodeCoords, double* physicalCoords) override;
        void getPhysicalCoords(const double* refCoords, const double* nodeCoords, double* physicalCoords) override;

    private:
        RefElem<ElemType::Triangle, BFType::Lagrangian,1> * geoElem; // Auxiliar reference element needed for geometry
    };


}

#endif //UPA_REFELEM_NEDELEC_H
