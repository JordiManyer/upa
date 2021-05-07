
#ifndef UPA_REFELEM_LAGRANGE_H
#define UPA_REFELEM_LAGRANGE_H

#include "elementDefinitions.h"
#include "quadratures.h"
#include "referenceElement.h"

namespace upa {

    /** @brief Templated class, derived from base class.
     *         Template specialisations will create different elements.
     */
    template<ElemType etype, int>
    class RefElem_Lagrange : public ReferenceElement {
    public:
        RefElem_Lagrange();
        ~RefElem_Lagrange() override = default;

        void evaluateBFs(const double* coords, double *bf) override;
        void evaluateDBFs(const double* coords, double *dbf) override;

        void getJacobian(int iG, const double* nodeCoords, double* J) override;
        void getJacobian(const double* dbf, const double* nodeCoords, double* J) override;

        void getPhysicalCoords(int iG, const double* nodeCoords, double* physicalCoords) override;
        void getPhysicalCoords(const double* refCoords, const double* nodeCoords, double* physicalCoords) override;

        void interpolateSolution(int iG, const double* dofs, double* sol) override;
        void interpolateSolution(const double* refCoords, const double* dofs, double* sol) override;

    };


    /*******************************************************************************************************************
     *****                                  TEMPLATE SPECIALISATIONS
     ******************************************************************************************************************/

    template <> RefElem_Lagrange<ElemType::Line,1>::RefElem_Lagrange();
    template <> void RefElem_Lagrange<ElemType::Line,1>::evaluateBFs(const double *coords, double *bf);
    template <> void RefElem_Lagrange<ElemType::Line,1>::evaluateDBFs(const double *coords, double *dbf);

    template <> RefElem_Lagrange<ElemType::Line,2>::RefElem_Lagrange();
    template <> void RefElem_Lagrange<ElemType::Line,2>::evaluateBFs(const double *coords, double *bf);
    template <> void RefElem_Lagrange<ElemType::Line,2>::evaluateDBFs(const double *coords, double *dbf);

    template <> RefElem_Lagrange<ElemType::Square,1>::RefElem_Lagrange();
    template <> void RefElem_Lagrange<ElemType::Square,1>::evaluateBFs(const double *coords, double *bf);
    template <> void RefElem_Lagrange<ElemType::Square,1>::evaluateDBFs(const double *coords, double *dbf);

    template <> RefElem_Lagrange<ElemType::Square,2>::RefElem_Lagrange();
    template <> void RefElem_Lagrange<ElemType::Square,2>::evaluateBFs(const double *coords, double *bf);
    template <> void RefElem_Lagrange<ElemType::Square,2>::evaluateDBFs(const double *coords, double *dbf);

    template <> RefElem_Lagrange<ElemType::Triangle,1>::RefElem_Lagrange();
    template <> void RefElem_Lagrange<ElemType::Triangle,1>::evaluateBFs(const double *coords, double *bf);
    template <> void RefElem_Lagrange<ElemType::Triangle,1>::evaluateDBFs(const double *coords, double *dbf);

    template <> RefElem_Lagrange<ElemType::Triangle,2>::RefElem_Lagrange();
    template <> void RefElem_Lagrange<ElemType::Triangle,2>::evaluateBFs(const double *coords, double *bf);
    template <> void RefElem_Lagrange<ElemType::Triangle,2>::evaluateDBFs(const double *coords, double *dbf);

}



#endif //UPA_REFELEM_LAGRANGE_H
