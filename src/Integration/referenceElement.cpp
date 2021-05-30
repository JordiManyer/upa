

#include "referenceElement.h"
#include "refElem_Lagrange.h"
#include "refElem_Nedelec.h"

#include <stdexcept>

namespace upa {


    ReferenceElement* getReferenceElement(ElemType etype, BFType bftype, int bforder) {
        if (etype == ElemType::Line) {
            if ( bftype == BFType::Lagrangian) {
                if (bforder == 1) return new RefElem_Lagrange<ElemType::Line,1>();
                if (bforder == 2) return new RefElem_Lagrange<ElemType::Line,2>();
            }
        } else if (etype == ElemType::Square) {
            if ( bftype == BFType::Lagrangian) {
                if (bforder == 1) return new RefElem_Lagrange<ElemType::Square,1>();
                if (bforder == 2) return new RefElem_Lagrange<ElemType::Square,2>();
            } else if ( bftype == BFType::Nedelec) {

            }
        } else if (etype == ElemType::Triangle) {
            if ( bftype == BFType::Lagrangian) {
                if (bforder == 1) return new RefElem_Lagrange<ElemType::Triangle,1>();
                if (bforder == 2) return new RefElem_Lagrange<ElemType::Triangle,2>();
            } else if ( bftype == BFType::Nedelec) {
                if (bforder == 1) return new RefElem_Nedelec<ElemType::Triangle,1>();
            }
        }
        throw std::runtime_error("getReferenceElement: Element not implemented!");
        return nullptr;
    }


    void ReferenceElement::evaluate(int degree, const double* refCoords, double *values) {
        if (degree == 0)      evaluateBFs(refCoords,values);
        else if (degree == 1) evaluateDBFs(refCoords,values);
        else throw std::runtime_error("ReferenceElement::evaluate : Derivative order too high!");
    }

}