
#include "refElem_Lagrangian.h"
#include "refElem_Nedelec.h"
#include <iostream>


namespace upa {


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


    RefElem<ElemType::Line, BFType::Lagrangian, 1>::RefElem() {
        _dim = 1;
        _nB = 2;
        _nSol = 1;
        _bforder = 1;
        _elemType = ElemType::Line;
        _bftype = BFType::Lagrangian;

        getQuadratureGauss_Line(1, _nG, &_gC, &_gW);

        double xi;
        _bf = new double[_nG * _nB];
        _dbf = new double[_nG * _nB * _dim];
        for (int iG = 0; iG < _nG; ++iG) {
            xi = _gC[iG * _dim + 0];

            _bf[iG * _nB + 0] = -(xi - 1) / 2;
            _bf[iG * _nB + 1] = (xi + 1) / 2;

            // dbf/dxi; dbf/deta
            _dbf[iG * _nB * _dim + 0 * _dim + 0] = -1.0 / 2;
            _dbf[iG * _nB * _dim + 1 * _dim + 0] = 1.0 / 2;

        }
    }

    void RefElem<ElemType::Line, BFType::Lagrangian, 1>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = -(refCoords[0]-1)/2;
        bf[1] = (refCoords[0]+1)/2;
    }

    void RefElem<ElemType::Line, BFType::Lagrangian, 1>::evaluateDBFs(const double *refCoords, double *dbf) {
        dbf[0] = -1.0/2.0;
        dbf[1] =  1.0/2.0;
    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


    RefElem<ElemType::Line, BFType::Lagrangian, 2>::RefElem() {
        _dim = 1;
        _nB = 3;
        _nSol = 1;
        _bforder = 2;
        _elemType = ElemType::Line;
        _bftype = BFType::Lagrangian;

        getQuadratureGauss_Line(2, _nG, &_gC, &_gW);

        double xi;
        _bf = new double[_nG * _nB];
        _dbf = new double[_nG * _nB * _dim];
        for (int iG = 0; iG < _nG; ++iG) {
            xi = _gC[iG * _dim + 0];

            _bf[iG * _nB + 0] = xi * (xi - 1.0) / 2.0;
            _bf[iG * _nB + 1] = 1.0 - xi * xi;
            _bf[iG * _nB + 2] = xi * (xi + 1.0) / 2.0;

            // dbf/dxi; dbf/deta
            _dbf[iG*_nB*_dim + 0*_dim + 0] = xi - 1.0 / 2.0;
            _dbf[iG*_nB*_dim + 1*_dim + 0] = - 2.0 * xi;
            _dbf[iG*_nB*_dim + 2*_dim + 0] = xi + 1.0 / 2.0;

        }
    }

    void RefElem<ElemType::Line, BFType::Lagrangian, 2>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = refCoords[0] * (refCoords[0] - 1.0) / 2.0;
        bf[1] = 1.0 - refCoords[0] * refCoords[0];
        bf[2] = refCoords[0] * (refCoords[0] + 1.0) / 2.0;
    }

    void RefElem<ElemType::Line, BFType::Lagrangian, 2>::evaluateDBFs(const double *refCoords, double *dbf) {
        dbf[0] =  refCoords[0] - 1.0 / 2.0;
        dbf[1] =  -2.0 * refCoords[0];
        dbf[2] =  refCoords[0] - 1.0 / 2.0;
    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


    RefElem<ElemType::Square, BFType::Lagrangian, 1>::RefElem() {
        _dim      = 2;
        _nB       = 4;
        _nSol     = 1;
        _bforder  = 1;
        _elemType = ElemType::Square;
        _bftype   = BFType::Lagrangian;

        getQuadratureGauss_Square(1, _nG, &_gC, &_gW);

        double xi,eta;
        _bf  = new double[_nG*_nB];
        _dbf = new double[_nG*_nB*_dim];
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

    void RefElem<ElemType::Square, BFType::Lagrangian, 1>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = (refCoords[0]-1)*(refCoords[1]-1)/4;
        bf[1] = -(refCoords[0]+1)*(refCoords[1]-1)/4;
        bf[2] = (refCoords[0]+1)*(refCoords[1]+1)/4;
        bf[3] = -(refCoords[0]-1)*(refCoords[1]+1)/4;
    }

    void RefElem<ElemType::Square, BFType::Lagrangian, 1>::evaluateDBFs(const double *refCoords, double *dbf) {
        dbf[0] = (refCoords[1]-1)/4  ;   dbf[1] = (refCoords[0]-1)/4 ;
        dbf[2] = -(refCoords[1]-1)/4 ;   dbf[3] = -(refCoords[0]+1)/4;
        dbf[4] = (refCoords[1]+1)/4  ;   dbf[5] = (refCoords[0]+1)/4 ;
        dbf[6] = -(refCoords[1]+1)/4 ;   dbf[7] = -(refCoords[0]-1)/4;
    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


    RefElem<ElemType::Square, BFType::Lagrangian, 2>::RefElem() {
        _dim      = 2;
        _nB       = 9;
        _nSol     = 1;
        _bforder  = 2;
        _elemType = ElemType::Square;
        _bftype   = BFType::Lagrangian;

        getQuadratureGauss_Square(2, _nG, &_gC, &_gW);

        double xi,eta;
        _bf  = new double[_nG*_nB];
        _dbf = new double[_nG*_nB*_dim];
        for (int iG = 0; iG < _nG; ++iG) {
            xi = _gC[iG*_dim+0]; eta = _gC[iG*_dim+1];

            _bf[iG*_nB + 0] = (xi*(xi-1.0)/2.0) * (eta*(eta-1.0)/2.0);
            _bf[iG*_nB + 1] = (xi*(xi-1.0)/2.0) * (1.0-eta*eta);
            _bf[iG*_nB + 2] = (xi*(xi-1.0)/2.0) * (eta*(eta+1.0)/2.0);
            
            _bf[iG*_nB + 3] = (1.0-xi*xi) * (eta*(eta-1.0)/2.0);
            _bf[iG*_nB + 4] = (1.0-xi*xi) * (1.0-eta*eta);
            _bf[iG*_nB + 5] = (1.0-xi*xi) * (eta*(eta+1.0)/2.0);
            
            _bf[iG*_nB + 6] = (xi*(xi+1.0)/2.0) * (eta*(eta-1.0)/2.0);
            _bf[iG*_nB + 7] = (xi*(xi+1.0)/2.0) * (1.0-eta*eta);
            _bf[iG*_nB + 8] = (xi*(xi+1.0)/2.0) * (eta*(eta+1.0)/2.0);


            // dbf/dxi; dbf/deta
            _dbf[iG*_nB*_dim + 0*_dim + 0] = (xi-1.0/2.0) * (eta*(eta-1.0)/2.0) ;  _dbf[iG*_nB*_dim + 0*_dim + 1] = (xi*(xi-1.0)/2.0) * (eta-1.0/2.0);
            _dbf[iG*_nB*_dim + 1*_dim + 0] = (xi-1.0/2.0) * (1.0-eta*eta);         _dbf[iG*_nB*_dim + 1*_dim + 1] = (xi*(xi-1.0)/2.0) * (- 2.0*eta);
            _dbf[iG*_nB*_dim + 2*_dim + 0] = (xi-1.0/2.0) * (eta*(eta+1.0)/2.0);   _dbf[iG*_nB*_dim + 2*_dim + 1] = (xi*(xi-1.0)/2.0) * (eta+1.0/2.0);

            _dbf[iG*_nB*_dim + 0*_dim + 0] = (-2.0*xi) * (eta*(eta-1.0)/2.0) ;  _dbf[iG*_nB*_dim + 0*_dim + 1] = (1.0-xi*xi) * (eta-1.0/2.0);
            _dbf[iG*_nB*_dim + 1*_dim + 0] = (-2.0*xi) * (1.0-eta*eta);         _dbf[iG*_nB*_dim + 1*_dim + 1] = (1.0-xi*xi) * (- 2.0*eta);
            _dbf[iG*_nB*_dim + 2*_dim + 0] = (-2.0*xi) * (eta*(eta+1.0)/2.0);   _dbf[iG*_nB*_dim + 2*_dim + 1] = (1.0-xi*xi) * (eta+1.0/2.0);

            _dbf[iG*_nB*_dim + 0*_dim + 0] = (xi+1.0/2.0) * (eta*(eta-1.0)/2.0) ;  _dbf[iG*_nB*_dim + 0*_dim + 1] = (xi*(xi+1.0)/2.0) * (eta-1.0/2.0);
            _dbf[iG*_nB*_dim + 1*_dim + 0] = (xi+1.0/2.0) * (1.0-eta*eta);         _dbf[iG*_nB*_dim + 1*_dim + 1] = (xi*(xi+1.0)/2.0) * (- 2.0*eta);
            _dbf[iG*_nB*_dim + 2*_dim + 0] = (xi+1.0/2.0) * (eta*(eta+1.0)/2.0);   _dbf[iG*_nB*_dim + 2*_dim + 1] = (xi*(xi+1.0)/2.0) * (eta+1.0/2.0);
            
        }
    }

    void RefElem<ElemType::Square, BFType::Lagrangian, 2>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = (refCoords[0]*(refCoords[0]-1.0)/2.0) * (refCoords[1]*(refCoords[1]-1.0)/2.0);
        bf[1] = (refCoords[0]*(refCoords[0]-1.0)/2.0) * (1.0-refCoords[1]*refCoords[1]);
        bf[2] = (refCoords[0]*(refCoords[0]-1.0)/2.0) * (refCoords[1]*(refCoords[1]+1.0)/2.0);

        bf[3] = (1.0-refCoords[0]*refCoords[0]) * (refCoords[1]*(refCoords[1]-1.0)/2.0);
        bf[4] = (1.0-refCoords[0]*refCoords[0]) * (1.0-refCoords[1]*refCoords[1]);
        bf[5] = (1.0-refCoords[0]*refCoords[0]) * (refCoords[1]*(refCoords[1]+1.0)/2.0);

        bf[6] = (refCoords[0]*(refCoords[0]+1.0)/2.0) * (refCoords[1]*(refCoords[1]-1.0)/2.0);
        bf[7] = (refCoords[0]*(refCoords[0]+1.0)/2.0) * (1.0-refCoords[1]*refCoords[1]);
        bf[8] = (refCoords[0]*(refCoords[0]+1.0)/2.0) * (refCoords[1]*(refCoords[1]+1.0)/2.0);
    }

    void RefElem<ElemType::Square, BFType::Lagrangian, 2>::evaluateDBFs(const double *refCoords, double *dbf) {
        // dN/dxi; dN/deta
        dbf[0] = (refCoords[0]-1.0/2.0) * (refCoords[1]*(refCoords[1]-1.0)/2.0) ;  dbf[1] = (refCoords[0]*(refCoords[0]-1.0)/2.0) * (refCoords[1]-1.0/2.0);
        dbf[2] = (refCoords[0]-1.0/2.0) * (1.0-refCoords[1]*refCoords[1]);         dbf[3] = (refCoords[0]*(refCoords[0]-1.0)/2.0) * (- 2.0*refCoords[1]);
        dbf[4] = (refCoords[0]-1.0/2.0) * (refCoords[1]*(refCoords[1]+1.0)/2.0);   dbf[5] = (refCoords[0]*(refCoords[0]-1.0)/2.0) * (refCoords[1]+1.0/2.0);

        dbf[6] = (-2.0*refCoords[0]) * (refCoords[1]*(refCoords[1]-1.0)/2.0) ;  dbf[7] = (1.0-refCoords[0]*refCoords[0]) * (refCoords[1]-1.0/2.0);
        dbf[8] = (-2.0*refCoords[0]) * (1.0-refCoords[1]*refCoords[1]);         dbf[9] = (1.0-refCoords[0]*refCoords[0]) * (- 2.0*refCoords[1]);
        dbf[10] = (-2.0*refCoords[0]) * (refCoords[1]*(refCoords[1]+1.0)/2.0);   dbf[11] = (1.0-refCoords[0]*refCoords[0]) * (refCoords[1]+1.0/2.0);

        dbf[12] = (refCoords[0]+1.0/2.0) * (refCoords[1]*(refCoords[1]-1.0)/2.0) ;  dbf[13] = (refCoords[0]*(refCoords[0]+1.0)/2.0) * (refCoords[1]-1.0/2.0);
        dbf[14] = (refCoords[0]+1.0/2.0) * (1.0-refCoords[1]*refCoords[1]);         dbf[15] = (refCoords[0]*(refCoords[0]+1.0)/2.0) * (- 2.0*refCoords[1]);
        dbf[16] = (refCoords[0]+1.0/2.0) * (refCoords[1]*(refCoords[1]+1.0)/2.0);   dbf[17] = (refCoords[0]*(refCoords[0]+1.0)/2.0) * (refCoords[1]+1.0/2.0);
    }

    ///**************************************************************************************************************///
    ///**************************************************************************************************************///

    /**  First order Lagrangian triangle:
     *
     *           N2                             N0 = 1-xi-eta
     *          [0,1]                           N1 = xi
     *            |   .                         N2 = eta
     *            |      .
     *            |         .
     *            |            .
     *          [0,0] -------- [1,0]
     *           N0              N1
     */

    RefElem<ElemType::Triangle, BFType::Lagrangian, 1>::RefElem() {
        _dim      = 2;
        _nB       = 3;
        _nSol     = 1;
        _bforder  = 1;
        _elemType = ElemType::Triangle;
        _bftype   = BFType::Lagrangian;

        getQuadratureGauss_Triangle(1, _nG, &_gC, &_gW);

        double xi,eta;
        _bf  = new double[_nG*_nB];
        _dbf = new double[_nG*_nB*_dim];
        for (int iG = 0; iG < _nG; ++iG) {
            xi = _gC[iG*_dim+0]; eta = _gC[iG*_dim+1];

            _bf[iG*_nB + 0] = 1.0 - xi -eta;
            _bf[iG*_nB + 1] = xi;
            _bf[iG*_nB + 2] = eta;

            // dbf/dxi; dbf/deta
            _dbf[iG*_nB*_dim + 0*_dim + 0] = -1.0;  _dbf[iG*_nB*_dim + 0*_dim + 1] = -1.0;
            _dbf[iG*_nB*_dim + 1*_dim + 0] = 1.0;   _dbf[iG*_nB*_dim + 1*_dim + 1] = 0.0;
            _dbf[iG*_nB*_dim + 2*_dim + 0] = 0.0;   _dbf[iG*_nB*_dim + 2*_dim + 1] = 1.0;

        }
    }

    void RefElem<ElemType::Triangle, BFType::Lagrangian, 1>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = 1.0 - refCoords[0] - refCoords[1];
        bf[1] = refCoords[0];
        bf[2] = refCoords[1];
    }

    void RefElem<ElemType::Triangle, BFType::Lagrangian, 1>::evaluateDBFs(const double *refCoords, double *dbf) {
        dbf[0] = -1.0  ; dbf[1] = -1.0 ;
        dbf[2] =  1.0  ; dbf[3] =  0.0 ;
        dbf[4] =  0.0  ; dbf[5] =  1.0 ;
    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///

    /**  Second order Lagrangian triangle:
     *
     *           N2                                  N0 = (1-xi-eta)(1-2 xi-2 eta)
     *          [0,1]                                N1 = xi (2 xi-1)
     *            |    .                             N2 = eta (2 eta-1)
     *            |        .                         N3 = 4 xi (1-xi-eta)
     *     N5  [0,1/2]     [1/2,1/2]   N4            N4 = 4 xi eta
     *            |              .                   N5 = 4 xi (1-xi-eta)
     *            |                  .
     *          [0,0] ---- [1/2,0] ---- [1,0]
     *           N0          N3           N1
     */

    RefElem<ElemType::Triangle, BFType::Lagrangian, 2>::RefElem() {
        _dim      = 2;
        _nB       = 9;
        _nSol     = 1;
        _bforder  = 2;
        _elemType = ElemType::Triangle;
        _bftype   = BFType::Lagrangian;

        getQuadratureGauss_Triangle(2, _nG, &_gC, &_gW);

        double xi,eta;
        _bf  = new double[_nG*_nB];
        _dbf = new double[_nG*_nB*_dim];
        for (int iG = 0; iG < _nG; ++iG) {
            xi = _gC[iG*_dim+0]; eta = _gC[iG*_dim+1];

            _bf[iG*_nB + 0] = (1.0-xi-eta) * (1-2.0*xi-2.0*eta);
            _bf[iG*_nB + 1] = xi*(2.0*xi-1.0);
            _bf[iG*_nB + 2] = eta*(2.0*xi-1.0);
            _bf[iG*_nB + 3] = 4.0*xi*(1.0-xi-eta);
            _bf[iG*_nB + 4] = 4.0*xi*eta;
            _bf[iG*_nB + 5] = 4.0*xi*(1.0-xi-eta);
            
            // dbf/dxi; dbf/deta
            _dbf[iG*_nB*_dim + 0*_dim + 0] = -3.0 + 4.0*xi + 4.0*eta;  _dbf[iG*_nB*_dim + 0*_dim + 1] = -3.0 + 4.0*xi + 4.0*eta;
            _dbf[iG*_nB*_dim + 1*_dim + 0] = 4.0*xi - 1.0;             _dbf[iG*_nB*_dim + 1*_dim + 1] = 0.0;
            _dbf[iG*_nB*_dim + 2*_dim + 0] = 2.0*eta;                  _dbf[iG*_nB*_dim + 2*_dim + 1] = 2.0*xi - 1.0;

            _dbf[iG*_nB*_dim + 0*_dim + 0] = 4.0 - 8.0*xi - 4.0*eta ;  _dbf[iG*_nB*_dim + 0*_dim + 1] = -4.0*xi;
            _dbf[iG*_nB*_dim + 1*_dim + 0] = 4.0*eta;                  _dbf[iG*_nB*_dim + 1*_dim + 1] = 4.0*xi;
            _dbf[iG*_nB*_dim + 2*_dim + 0] = 4.0 - 8.0*xi - 4.0*eta;   _dbf[iG*_nB*_dim + 2*_dim + 1] = -4.0*xi;
            
        }
    }

    void RefElem<ElemType::Triangle, BFType::Lagrangian, 2>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = (1.0-refCoords[0]-refCoords[1]) * (1-2.0*refCoords[0]-2.0*refCoords[1]);
        bf[1] = refCoords[0]*(2.0*refCoords[0]-1.0);
        bf[2] = refCoords[1]*(2.0*refCoords[0]-1.0);
        bf[3] = 4.0*refCoords[0]*(1.0-refCoords[0]-refCoords[1]);
        bf[4] = 4.0*refCoords[0]*refCoords[1];
        bf[5] = 4.0*refCoords[0]*(1.0-refCoords[0]-refCoords[1]);
    }

    void RefElem<ElemType::Triangle, BFType::Lagrangian, 2>::evaluateDBFs(const double *refCoords, double *dbf) {
        // dN/dxi; dN/deta
        dbf[0] = -3.0 + 4.0*refCoords[0] + 4.0*refCoords[1];  dbf[1] = -3.0 + 4.0*refCoords[0] + 4.0*refCoords[1];
        dbf[2] = 4.0*refCoords[0] - 1.0;                   dbf[3] = 0.0;
        dbf[4] = 2.0*refCoords[1];                         dbf[5] = 2.0*refCoords[0] - 1.0;

        dbf[6]  = 4.0 - 8.0*refCoords[0] - 4.0*refCoords[1] ; dbf[7]  = -4.0*refCoords[0];
        dbf[8]  = 4.0*refCoords[1];                        dbf[9]  = 4.0*refCoords[0];
        dbf[10] = 4.0 - 8.0*refCoords[0] - 4.0*refCoords[1];  dbf[11] = -4.0*refCoords[0];
    }

    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


}