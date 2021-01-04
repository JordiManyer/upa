

#include "referenceElement.h"

namespace upa {

    ReferenceElement::ReferenceElement(ElemType etype, BFType bftype, int bforder) {
        if (bftype == BFType::Lagrangian) {

            if      (etype == ElemType::Line  ) _fillRefElem<ElemType::Line, BFType::Lagrangian>(bforder);
            else if (etype == ElemType::Square) _fillRefElem<ElemType::Square, BFType::Lagrangian>(bforder);
            else throw std::runtime_error("ReferenceElement: ElemType not implemented.");

        } else throw std::runtime_error("ReferenceElement: BFtype not implemented.");
    }



    template <>
    void ReferenceElement::_fillRefElem<ElemType::Line, BFType::Lagrangian>(int bforder) {

        _dim      = 1;
        _nB       = 2;
        _bforder  = bforder;
        _elemType = ElemType::Line;
        _bftype   = BFType::Lagrangian;

        switch (_bforder) {
            case 1:
                _nG  = 2;
                _gW  = new double[_nG];
                _gC  = new double[_nG*_dim];
                _bf  = new double[_nG*_nB];
                _dbf = new double[_nG*_nB*_dim];

                _gW[0] = 1.0; _gW[1] = 1.0;
                _gC[0] = -1/sqrt(3); _gC[1] = 1/sqrt(3);

                double xi;
                for (int iG = 0; iG < _nG; ++iG) {
                    xi = _gC[iG*_dim+0];

                    _bf[iG*_nB + 0] = -(xi-1)/2;
                    _bf[iG*_nB + 1] = (xi+1)/2;

                    // dbf/dxi; dbf/deta
                    _dbf[iG*_nB*_dim + 0*_dim + 0] = -1.0/2 ;
                    _dbf[iG*_nB*_dim + 1*_dim + 0] = 1.0/2;
                }
                break;

            default:
                throw std::runtime_error("ReferenceElement: BF order not implemented.");
        }
    }


    template <>
    void ReferenceElement::_fillRefElem<ElemType::Square, BFType::Lagrangian>(int bforder) {

        _dim      = 2;
        _nB       = 4;
        _bforder  = bforder;
        _elemType = ElemType::Square;
        _bftype   = BFType::Lagrangian;

        switch (_bforder) {
            case 1:
                _nG  = 4;
                _gW  = new double[_nG];
                _gC  = new double[_nG*_dim];
                _bf  = new double[_nG*_nB];
                _dbf = new double[_nG*_nB*_dim];

                _gW[0] = 1.0; _gW[1] = 1.0; _gW[2] = 1.0; _gW[3] = 1.0;
                _gC[0] = -1/sqrt(3); _gC[2] = 1/sqrt(3); _gC[4] = 1/sqrt(3); _gC[6] = -1/sqrt(3); // x
                _gC[1] = -1/sqrt(3); _gC[3] = -1/sqrt(3); _gC[5] = 1/sqrt(3); _gC[7] = 1/sqrt(3); // y

                double xi,eta;
                for (int iG = 0; iG < _nG; ++iG) {
                    xi = _gC[iG*_dim+0]; eta = _gC[iG*_dim+1];

                    _bf[iG*_nB + 0] = (xi-1)*(eta-1)/4;
                    _bf[iG*_nB + 1] = -(xi+1)*(eta-1)/4;
                    _bf[iG*_nB + 2] = -(xi-1)*(eta+1)/4;
                    _bf[iG*_nB + 3] = (xi+1)*(eta+1)/4;

                    // dbf/dxi; dbf/deta
                    _dbf[iG*_nB*_dim + 0*_dim + 0] = (eta-1)/4 ;   _dbf[iG*_nB*_dim + 0*_dim + 1] = (xi-1)/4 ;
                    _dbf[iG*_nB*_dim + 1*_dim + 0] = -(eta-1)/4;   _dbf[iG*_nB*_dim + 1*_dim + 1] = -(xi+1)/4;
                    _dbf[iG*_nB*_dim + 2*_dim + 0] = -(eta+1)/4;   _dbf[iG*_nB*_dim + 2*_dim + 1] = -(xi-1)/4;
                    _dbf[iG*_nB*_dim + 3*_dim + 0] = (eta+1)/4 ;   _dbf[iG*_nB*_dim + 3*_dim + 1] = (xi+1)/4 ;
                }
                break;

            default:
                throw std::runtime_error("ReferenceElement: BF order not implemented.");
        }
    }


}