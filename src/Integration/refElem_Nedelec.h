
#ifndef UPA_REFELEM_NEDELEC_H
#define UPA_REFELEM_NEDELEC_H

#include "elementDefinitions.h"
#include "referenceElement.h"
#include "refElem_Lagrange.h"

namespace upa {

    /** @brief Templated class, derived from base class.
     *         Template specialisations will create different elements.
     */
    template<ElemType etype, int>
    class RefElem_Nedelec : public ReferenceElement {
    public:
        RefElem_Nedelec();
        ~RefElem_Nedelec() override = default;

        void evaluateBFs(const double* coords, double *bf) override;
        void evaluateDBFs(const double* coords, double *dbf) override;

        void getJacobian(int iG, const double* nodeCoords, double* J) override;
        void getJacobian(const double* dbf, const double* nodeCoords, double* J) override;

        void getPhysicalCoords(int iG, const double* nodeCoords, double* physicalCoords) override;
        void getPhysicalCoords(const double* refCoords, const double* nodeCoords, double* physicalCoords) override;

        void interpolateSolution(int iG, const double* dofs, double* sol) override;
        void interpolateSolution(const double* refCoords, const double* dofs, double* sol) override;

    private:
        // Lagrangian basis functions for geometry
        double *_geo_bf;
        double *_geo_dbf;

        void geo_evaluateBFs(const double* refCoords, double *bf);
        void geo_evaluateDBFs(const double* refCoords, double *dbf);
    };



    /*******************************************************************************************************************
     *****                                  TEMPLATE SPECIALISATIONS
     ******************************************************************************************************************/

    template <> RefElem_Nedelec<ElemType::Triangle, 1>::RefElem_Nedelec();
    template <> void RefElem_Nedelec<ElemType::Triangle, 1>::evaluateBFs(const double *coords, double *bf);
    template <> void RefElem_Nedelec<ElemType::Triangle, 1>::evaluateDBFs(const double *coords, double *dbf);
    template <> void RefElem_Nedelec<ElemType::Triangle, 1>::geo_evaluateBFs(const double *refCoords, double *bf);
    template <> void RefElem_Nedelec<ElemType::Triangle, 1>::geo_evaluateDBFs(const double *refCoords, double *dbf);

}

#endif //UPA_REFELEM_NEDELEC_H
