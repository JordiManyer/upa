
#include "refElem_Nedelec.h"

namespace upa {

    /*******************************************************************************************************************
     *****                                    GENERAL FUNCTIONS
     ******************************************************************************************************************/
    template<ElemType etype, int order>
    RefElem_Nedelec<etype, order>::RefElem_Nedelec() {
        throw std::runtime_error("RefElem_Lagrange: Element not implemented!");
    }

    template<ElemType etype, int order>
    void RefElem_Nedelec<etype, order>::getJacobian(int iG, const double *nodeCoords, double *J) {
        for (int i = 0; i < _dim; ++i)
            for (int j = 0; j < _dim; ++j) {
                J[i*_dim + j] = 0.0;
                for (int l = 0; l < _nB; ++l)
                    J[i*_dim + j] += _geo_dbf[_nB*_dim*iG + _dim*l + i] * nodeCoords[_dim*l + j];
            }
    }

    template<ElemType etype, int order>
    void RefElem_Nedelec<etype, order>::getJacobian(const double* refCoords, const double* nodeCoords, double* J) {
        double dbf[_nB*_dim];
        geo_evaluateDBFs(refCoords,dbf);
        for (int i = 0; i < _dim; ++i)
            for (int j = 0; j < _dim; ++j) {
                J[i*_dim + j] = 0.0;
                for (int l = 0; l < _nB; ++l)
                    J[i*_dim + j] += dbf[_dim*l + i] * nodeCoords[_dim*l + j];
            }
    }

    template<ElemType etype, int order>
    void RefElem_Nedelec<etype, order>::getPhysicalCoords(int iG, const double* nodeCoords, double* physicalCoords) {
        for (int i= 0; i < _dim; ++i) {
            physicalCoords[i] = 0.0;
            for (int j = 0; j < _nB; ++j) physicalCoords[i] += nodeCoords[j*_dim + i] * _geo_bf[iG*_nB + j];
        }
    }

    template<ElemType etype, int order>
    void RefElem_Nedelec<etype, order>::getPhysicalCoords(const double* refCoords, const double* nodeCoords, double* physicalCoords) {
        double bf[_nB];
        geo_evaluateBFs(refCoords,bf);
        for (int i= 0; i < _dim; ++i) {
            physicalCoords[i] = 0.0;
            for (int j = 0; j < _nB; ++j) physicalCoords[i] += nodeCoords[j*_dim + i] * bf[j];
        }
    }

    template<ElemType etype, int order>
    void RefElem_Nedelec<etype, order>::interpolateSolution(int iG, const double* dofs, double* sol) {
        for (int i = 0; i < _nSol; ++i) {
            sol[i] = 0.0;
            for (int j = 0; j < _nB; ++j) sol[i] += _bf[iG*_nB*_nSol + j*_nSol + i] * dofs[j];
        }
    }

    template<ElemType etype, int order>
    void RefElem_Nedelec<etype, order>::interpolateSolution(const double* refCoords, const double* dofs, double* sol) {
        double bf[_nB * _nSol];
        evaluateBFs(refCoords, bf);
        for (int i = 0; i < _nSol; ++i) {
            sol[i] = 0.0;
            for (int j = 0; j < _nB; ++j) sol[i] += bf[j *_nSol + i] * dofs[j];
        }
    }

    ///**************************************************************************************************************///
    ///**************************************************************************************************************///

    /**  First order Lagrangian triangle:
     *   We have one DOF associated to each edge ei, with a basis function Ni such that
     *                   integral_ei { ti·Nj ds } = delta_ij,
     *   with ti the unit tangent to the edge ei.
     *
     *                                          N0 = [1-eta , xi]
     *          [0,1]                           N1 = [-eta , xi]
     *            |   .                         N2 = [-eta , xi - 1]
     *            |      .  e1
     *        e2  |         .
     *            |            .
     *          [0,0] -------- [1,0]
     *                    e0
     */

    template <>
    RefElem_Nedelec<ElemType::Triangle, 1>::RefElem_Nedelec() {
        _dim = 2;
        _nB = 3;
        _nSol = 2;
        _bforder = 1;
        _elemType = ElemType::Triangle;
        _bftype = BFType::Nedelec;

        getQuadratureGauss_Triangle(1, _nG, &_gC, &_gW);

        double xi, eta;
        _bf      = new double[_nG * _nB * _dim];
        _dbf     = new double[_nG * _nB * _dim * _dim];
        _curlbf  = new double[_nG * _nB];
        _geo_bf  = new double[_nG * _nB];
        _geo_dbf = new double[_nG * _nB * _dim];
        for (int iG = 0; iG < _nG; ++iG) {
            xi = _gC[iG * _dim + 0];
            eta = _gC[iG * _dim + 1];

            _bf[iG*_nB*_nSol + 0*_nSol + 0] = 1.0 - eta; _bf[iG*_nB*_nSol + 0*_nSol + 1] = xi;
            _bf[iG*_nB*_nSol + 1*_nSol + 0] = - eta;     _bf[iG*_nB*_nSol + 1*_nSol + 1] = xi;
            _bf[iG*_nB*_nSol + 2*_nSol + 0] = - eta;     _bf[iG*_nB*_nSol + 2*_nSol + 1] = xi - 1.0;

            // dbf/dxi; dbf/deta
            _dbf[iG*_nB*_nSol*_dim + 0*_nSol*_dim + 0*_dim + 0] = 0.0; _dbf[iG*_nB*_nSol*_dim + 0*_nSol*_dim + 0*_dim + 1] = -1.0;
            _dbf[iG*_nB*_nSol*_dim + 0*_nSol*_dim + 1*_dim + 0] = 1.0; _dbf[iG*_nB*_nSol*_dim + 0*_nSol*_dim + 1*_dim + 1] = 0.0;

            _dbf[iG*_nB*_nSol*_dim + 1*_nSol*_dim + 0*_dim + 0] = 0.0; _dbf[iG*_nB*_nSol*_dim + 1*_nSol*_dim + 0*_dim + 1] = -1.0;
            _dbf[iG*_nB*_nSol*_dim + 1*_nSol*_dim + 1*_dim + 0] = 1.0; _dbf[iG*_nB*_nSol*_dim + 1*_nSol*_dim + 1*_dim + 1] = 0.0;

            _dbf[iG*_nB*_nSol*_dim + 2*_nSol*_dim + 0*_dim + 0] = 0.0; _dbf[iG*_nB*_nSol*_dim + 2*_nSol*_dim + 0*_dim + 1] = -1.0;
            _dbf[iG*_nB*_nSol*_dim + 2*_nSol*_dim + 1*_dim + 0] = 1.0; _dbf[iG*_nB*_nSol*_dim + 2*_nSol*_dim + 1*_dim + 1] = 0.0;

            _curlbf[iG*_nB + 0] = 2.0;
            _curlbf[iG*_nB + 1] = 2.0;
            _curlbf[iG*_nB + 2] = 2.0;

            // Geometry
            _geo_bf[iG*_nB + 0] = 1.0 - xi -eta;
            _geo_bf[iG*_nB + 1] = xi;
            _geo_bf[iG*_nB + 2] = eta;

            _geo_dbf[iG*_nB*_dim + 0*_dim + 0] = -1.0;  _geo_dbf[iG*_nB*_dim + 0*_dim + 1] = -1.0;
            _geo_dbf[iG*_nB*_dim + 1*_dim + 0] = 1.0;   _geo_dbf[iG*_nB*_dim + 1*_dim + 1] = 0.0;
            _geo_dbf[iG*_nB*_dim + 2*_dim + 0] = 0.0;   _geo_dbf[iG*_nB*_dim + 2*_dim + 1] = 1.0;
        }
    }

    template <>
    void RefElem_Nedelec<ElemType::Triangle, 1>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = 1.0 - refCoords[1]; bf[1] = refCoords[0];
        bf[2] = - refCoords[1];     bf[3] = refCoords[0];
        bf[4] = - refCoords[1];     bf[5] = refCoords[0] - 1.0;
    }

    template <>
    void RefElem_Nedelec<ElemType::Triangle, 1>::evaluateDBFs(const double *refCoords, double *dbf) {
        dbf[0] = 0.0; dbf[1] = -1.0;
        dbf[2] = 1.0; dbf[3] = 0.0;
        dbf[4] = 0.0; dbf[5] = -1.0;
        dbf[6] = 1.0; dbf[7] = 0.0;
        dbf[8] = 0.0; dbf[9] = -1.0;
        dbf[10] = 1.0;dbf[11] = 0.0;
    }

    template <>
    void RefElem_Nedelec<ElemType::Triangle, 1>::geo_evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = 1.0 - refCoords[0] - refCoords[1];
        bf[1] = refCoords[0];
        bf[2] = refCoords[1];
    }

    template <>
    void RefElem_Nedelec<ElemType::Triangle, 1>::geo_evaluateDBFs(const double *refCoords, double *dbf) {
        dbf[0] = -1.0  ; dbf[1] = -1.0 ;
        dbf[2] =  1.0  ; dbf[3] =  0.0 ;
        dbf[4] =  0.0  ; dbf[5] =  1.0 ;
    }

    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


}