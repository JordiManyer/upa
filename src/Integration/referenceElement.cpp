

#include "referenceElement.h"
#include "refElem_Lagrangian.h"
#include "refElem_Nedelec.h"

namespace upa {


    ReferenceElement* getReferenceElement(ElemType etype, BFType bftype, int bforder) {
        if (etype == ElemType::Line) {
            if ( bftype == BFType::Lagrangian) {
                if (bforder == 1) return new RefElem<ElemType::Line,BFType::Lagrangian,1>();
                if (bforder == 2) return new RefElem<ElemType::Line,BFType::Lagrangian,2>();
            } else if ( bftype == BFType::Nedelec) {

            }
        } else if (etype == ElemType::Square) {
            if ( bftype == BFType::Lagrangian) {
                if (bforder == 1) return new RefElem<ElemType::Square,BFType::Lagrangian,1>();
                if (bforder == 2) return new RefElem<ElemType::Square,BFType::Lagrangian,2>();
            } else if ( bftype == BFType::Nedelec) {

            }
        } else if (etype == ElemType::Triangle) {
            if ( bftype == BFType::Lagrangian) {
                if (bforder == 1) return new RefElem<ElemType::Triangle,BFType::Lagrangian,1>();
                if (bforder == 2) return new RefElem<ElemType::Triangle,BFType::Lagrangian,2>();
            } else if ( bftype == BFType::Nedelec) {
                if (bforder == 1) return new RefElem<ElemType::Triangle,BFType::Nedelec,1>();
            }
        }
        throw std::runtime_error("getReferenceElement: Element not implemented!");
        return nullptr;
    }


    void ReferenceElement::evaluate(int degree, const double* coords, double *values) {
        if (degree == 0)      evaluateBFs(coords,values);
        else if (degree == 1) evaluateDBFs(coords,values);
        else throw std::runtime_error("ReferenceElement::evaluate : Derivative order too high!");
    }


}