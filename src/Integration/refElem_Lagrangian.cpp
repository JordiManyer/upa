
#include "refElem_Lagrangian.h"



namespace upa {


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


    RefElem<ElemType::Line, BFType::Lagrangian, 1>::RefElem() {
        _dim = 1;
        _nB = 2;
        _bforder = 1;
        _elemType = ElemType::Line;
        _bftype = BFType::Lagrangian;

        _nG = 2;
        _gW = new double[_nG];
        _gC = new double[_nG * _dim];
        _bf = new double[_nG * _nB];
        _dbf = new double[_nG * _nB * _dim];

        _gW[0] = 1.0; _gW[1] = 1.0;
        _gC[0] = -1 / sqrt(3); _gC[1] = 1 / sqrt(3);

        double xi;
        for (int iG = 0; iG < _nG; ++iG) {
            xi = _gC[iG * _dim + 0];

            _bf[iG * _nB + 0] = -(xi - 1) / 2;
            _bf[iG * _nB + 1] = (xi + 1) / 2;

            // dbf/dxi; dbf/deta
            _dbf[iG * _nB * _dim + 0 * _dim + 0] = -1.0 / 2;
            _dbf[iG * _nB * _dim + 1 * _dim + 0] = 1.0 / 2;

        }
    }

    void RefElem<ElemType::Line, BFType::Lagrangian, 1>::evaluateBFs(const double *coords, double *bf) {
        bf[0] = -(coords[0]-1)/2;
        bf[1] = (coords[0]+1)/2;
    }

    void RefElem<ElemType::Line, BFType::Lagrangian, 1>::evaluateDBFs(const double *coords, double *dbf) {
        dbf[0] = -1.0/2.0;
        dbf[1] =  1.0/2.0;
    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


    RefElem<ElemType::Line, BFType::Lagrangian, 2>::RefElem() {

    }

    void RefElem<ElemType::Line, BFType::Lagrangian, 2>::evaluateBFs(const double *coords, double *bf) {

    }

    void RefElem<ElemType::Line, BFType::Lagrangian, 2>::evaluateDBFs(const double *coords, double *dbf) {

    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


    RefElem<ElemType::Square, BFType::Lagrangian, 1>::RefElem() {
        _dim      = 2;
        _nB       = 4;
        _bforder  = 1;
        _elemType = ElemType::Square;
        _bftype   = BFType::Lagrangian;

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
            _bf[iG*_nB + 2] = (xi+1)*(eta+1)/4;
            _bf[iG*_nB + 3] = -(xi-1)*(eta+1)/4;

            // dbf/dxi; dbf/deta
            _dbf[iG*_nB*_dim + 0*_dim + 0] = (eta-1)/4 ;   _dbf[iG*_nB*_dim + 0*_dim + 1] = (xi-1)/4 ;
            _dbf[iG*_nB*_dim + 1*_dim + 0] = -(eta-1)/4;   _dbf[iG*_nB*_dim + 1*_dim + 1] = -(xi+1)/4;
            _dbf[iG*_nB*_dim + 2*_dim + 0] = (eta+1)/4;   _dbf[iG*_nB*_dim + 2*_dim + 1] = (xi+1)/4;
            _dbf[iG*_nB*_dim + 3*_dim + 0] = -(eta+1)/4 ;   _dbf[iG*_nB*_dim + 3*_dim + 1] = -(xi-1)/4 ;
        }
    }

    void RefElem<ElemType::Square, BFType::Lagrangian, 1>::evaluateBFs(const double *coords, double *bf) {
        bf[0] = (coords[0]-1)*(coords[1]-1)/4;
        bf[1] = -(coords[0]+1)*(coords[1]-1)/4;
        bf[2] = (coords[0]+1)*(coords[1]+1)/4;
        bf[3] = -(coords[0]-1)*(coords[1]+1)/4;
    }

    void RefElem<ElemType::Square, BFType::Lagrangian, 1>::evaluateDBFs(const double *coords, double *dbf) {
        dbf[0] = (coords[1]-1)/4  ;   dbf[1] = (coords[0]-1)/4 ;
        dbf[2] = -(coords[1]-1)/4 ;   dbf[3] = -(coords[0]+1)/4;
        dbf[4] = (coords[1]+1)/4  ;   dbf[5] = (coords[0]+1)/4 ;
        dbf[6] = -(coords[1]+1)/4 ;   dbf[7] = -(coords[0]-1)/4;
    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


    RefElem<ElemType::Square, BFType::Lagrangian, 2>::RefElem() {

    }

    void RefElem<ElemType::Square, BFType::Lagrangian, 2>::evaluateBFs(const double *coords, double *bf) {

    }

    void RefElem<ElemType::Square, BFType::Lagrangian, 2>::evaluateDBFs(const double *coords, double *dbf) {

    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


}