
#include "refElem_Lagrange.h"
#include "quadratures.h"

namespace upa {

    /*******************************************************************************************************************
     *****                                    GENERAL FUNCTIONS
     ******************************************************************************************************************/
    template<ElemType etype, int order>
    RefElem_Lagrange<etype, order>::RefElem_Lagrange() {
        throw std::runtime_error("RefElem_Lagrange: Element not implemented!");
    }

    template<ElemType etype, int order>
    void RefElem_Lagrange<etype, order>::getJacobian(int iG, const double *nodeCoords, double *J) {
        for (int i = 0; i < _dim; ++i)
            for (int j = 0; j < _dim; ++j) {
                J[i*_dim + j] = 0.0;
                for (int l = 0; l < _nB; ++l)
                    J[i*_dim + j] += _dbf[_nB*_dim*iG + _dim*l + i] * nodeCoords[_dim*l + j];
            }
    }

    template<ElemType etype, int order>
    void RefElem_Lagrange<etype, order>::getJacobian(const double* refCoords, const double* nodeCoords, double* J) {
        double dbf[_nB*_dim];
        evaluateDBFs(refCoords,dbf);
        for (int i = 0; i < _dim; ++i)
            for (int j = 0; j < _dim; ++j) {
                J[i*_dim + j] = 0.0;
                for (int l = 0; l < _nB; ++l)
                    J[i*_dim + j] += dbf[_dim*l + i] * nodeCoords[_dim*l + j];
            }
    }

    template<ElemType etype, int order>
    void RefElem_Lagrange<etype, order>::getPhysicalCoords(int iG, const double* nodeCoords, double* physicalCoords) {
        for (int i= 0; i < _dim; ++i) {
            physicalCoords[i] = 0.0;
            for (int j = 0; j < _nB; ++j) physicalCoords[i] += nodeCoords[j*_dim + i] * _bf[iG*_nB + j];
        }
    }

    template<ElemType etype, int order>
    void RefElem_Lagrange<etype, order>::getPhysicalCoords(const double* refCoords, const double* nodeCoords, double* physicalCoords) {
        double bf[_nB];
        evaluateBFs(refCoords,bf);
        for (int i= 0; i < _dim; ++i) {
            physicalCoords[i] = 0.0;
            for (int j = 0; j < _nB; ++j) physicalCoords[i] += nodeCoords[j*_dim + i] * bf[j];
        }
    }

    template<ElemType etype, int order>
    void RefElem_Lagrange<etype, order>::interpolateSolution(int iG, const double* dofs, double* sol) {
        *sol = 0.0;
        for (int j = 0; j < _nB; ++j) *sol += _bf[iG*_nB + j] * dofs[j];
    }

    template<ElemType etype, int order>
    void RefElem_Lagrange<etype, order>::interpolateSolution(const double* refCoords, const double* dofs, double* sol) {
        double bf[_nB];
        evaluateBFs(refCoords, bf);
        *sol = 0.0;
        for (int j = 0; j < _nB; ++j) *sol += bf[j] * dofs[j];
    }

    
    
    ///**************************************************************************************************************///
    ///                                      SPECIALISATIONS                                                         ///
    ///**************************************************************************************************************///

    template <>
    RefElem_Lagrange<ElemType::Line, 1>::RefElem_Lagrange() {
        _dim = 1;
        _nB = 2;
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

    template <>
    void RefElem_Lagrange<ElemType::Line, 1>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = -(refCoords[0]-1)/2;
        bf[1] = (refCoords[0]+1)/2;
    }

    template <>
    void RefElem_Lagrange<ElemType::Line, 1>::evaluateDBFs(const double *refCoords, double *dbf) {
        dbf[0] = -1.0/2.0;
        dbf[1] =  1.0/2.0;
    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///

    template <>
    RefElem_Lagrange<ElemType::Line, 2>::RefElem_Lagrange() {
        _dim = 1;
        _nB = 3;
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

    template <>
    void RefElem_Lagrange<ElemType::Line, 2>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = refCoords[0] * (refCoords[0] - 1.0) / 2.0;
        bf[1] = 1.0 - refCoords[0] * refCoords[0];
        bf[2] = refCoords[0] * (refCoords[0] + 1.0) / 2.0;
    }

    template<>
    void RefElem_Lagrange<ElemType::Line, 2>::evaluateDBFs(const double *refCoords, double *dbf) {
        dbf[0] =  refCoords[0] - 1.0 / 2.0;
        dbf[1] =  -2.0 * refCoords[0];
        dbf[2] =  refCoords[0] - 1.0 / 2.0;
    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///

    template <>
    RefElem_Lagrange<ElemType::Square, 1>::RefElem_Lagrange() {
        _dim      = 2;
        _nB       = 4;
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

    template <>
    void RefElem_Lagrange<ElemType::Square, 1>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = (refCoords[0]-1)*(refCoords[1]-1)/4;
        bf[1] = -(refCoords[0]+1)*(refCoords[1]-1)/4;
        bf[2] = (refCoords[0]+1)*(refCoords[1]+1)/4;
        bf[3] = -(refCoords[0]-1)*(refCoords[1]+1)/4;
    }

    template <>
    void RefElem_Lagrange<ElemType::Square, 1>::evaluateDBFs(const double *refCoords, double *dbf) {
        dbf[0] = (refCoords[1]-1)/4  ;   dbf[1] = (refCoords[0]-1)/4 ;
        dbf[2] = -(refCoords[1]-1)/4 ;   dbf[3] = -(refCoords[0]+1)/4;
        dbf[4] = (refCoords[1]+1)/4  ;   dbf[5] = (refCoords[0]+1)/4 ;
        dbf[6] = -(refCoords[1]+1)/4 ;   dbf[7] = -(refCoords[0]-1)/4;
    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///

    template <>
    RefElem_Lagrange<ElemType::Square, 2>::RefElem_Lagrange() {
        _dim      = 2;
        _nB       = 9;
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

    template <>
    void RefElem_Lagrange<ElemType::Square, 2>::evaluateBFs(const double *refCoords, double *bf) {
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

    template <>
    void RefElem_Lagrange<ElemType::Square, 2>::evaluateDBFs(const double *refCoords, double *dbf) {
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

    template <>
    RefElem_Lagrange<ElemType::Triangle, 1>::RefElem_Lagrange() {
        _dim      = 2;
        _nB       = 3;
        _bforder  = 1;
        _elemType = ElemType::Triangle;
        _bftype   = BFType::Lagrangian;

        getQuadratureGauss_Triangle(2, _nG, &_gC, &_gW);

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

    template <>
    void RefElem_Lagrange<ElemType::Triangle, 1>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = 1.0 - refCoords[0] - refCoords[1];
        bf[1] = refCoords[0];
        bf[2] = refCoords[1];
    }

    template <>
    void RefElem_Lagrange<ElemType::Triangle, 1>::evaluateDBFs(const double *refCoords, double *dbf) {
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

    template <>
    RefElem_Lagrange<ElemType::Triangle, 2>::RefElem_Lagrange() {
        _dim      = 2;
        _nB       = 9;
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

    template <>
    void RefElem_Lagrange<ElemType::Triangle, 2>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = (1.0-refCoords[0]-refCoords[1]) * (1-2.0*refCoords[0]-2.0*refCoords[1]);
        bf[1] = refCoords[0]*(2.0*refCoords[0]-1.0);
        bf[2] = refCoords[1]*(2.0*refCoords[0]-1.0);
        bf[3] = 4.0*refCoords[0]*(1.0-refCoords[0]-refCoords[1]);
        bf[4] = 4.0*refCoords[0]*refCoords[1];
        bf[5] = 4.0*refCoords[0]*(1.0-refCoords[0]-refCoords[1]);
    }

    template <>
    void RefElem_Lagrange<ElemType::Triangle, 2>::evaluateDBFs(const double *refCoords, double *dbf) {
        // dN/dxi; dN/deta
        dbf[0] = -3.0 + 4.0*refCoords[0] + 4.0*refCoords[1];  dbf[1] = -3.0 + 4.0*refCoords[0] + 4.0*refCoords[1];
        dbf[2] = 4.0*refCoords[0] - 1.0;                      dbf[3] = 0.0;
        dbf[4] = 2.0*refCoords[1];                            dbf[5] = 2.0*refCoords[0] - 1.0;

        dbf[6]  = 4.0 - 8.0*refCoords[0] - 4.0*refCoords[1] ; dbf[7]  = -4.0*refCoords[0];
        dbf[8]  = 4.0*refCoords[1];                           dbf[9]  = 4.0*refCoords[0];
        dbf[10] = 4.0 - 8.0*refCoords[0] - 4.0*refCoords[1];  dbf[11] = -4.0*refCoords[0];
    }

    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


}