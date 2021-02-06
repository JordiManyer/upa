

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


    void ReferenceElement::evaluate(int degree, const double* refCoords, double *values) {
        if (degree == 0)      evaluateBFs(refCoords,values);
        else if (degree == 1) evaluateDBFs(refCoords,values);
        else throw std::runtime_error("ReferenceElement::evaluate : Derivative order too high!");
    }


    void ReferenceElement::getJacobian(int iG, const double *nodeCoords, double *J) {
        for (int i = 0; i < _dim; ++i)
            for (int j = 0; j < _dim; ++j) {
                J[i*_dim + j] = 0.0;
                for (int l = 0; l < _nB; ++l)
                    J[i*_dim + j] += _dbf[_nB*_dim*iG + _dim*l + i] * nodeCoords[_dim*l + j];
            }
    }


    void ReferenceElement::getJacobian(const double* refCoords, const double* nodeCoords, double* J) {
        double dbf[_nB*_dim];
        evaluateDBFs(refCoords,dbf);
        for (int i = 0; i < _dim; ++i)
            for (int j = 0; j < _dim; ++j) {
                J[i*_dim + j] = 0.0;
                for (int l = 0; l < _nB; ++l)
                    J[i*_dim + j] += dbf[_dim*l + i] * nodeCoords[_dim*l + j];
            }
    }

}