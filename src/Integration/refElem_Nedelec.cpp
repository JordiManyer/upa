
#include "refElem_Nedelec.h"

namespace upa {

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

    RefElem<ElemType::Triangle, BFType::Nedelec, 1>::RefElem() {
        _dim = 2;
        _nB = 3;
        _bforder = 1;
        _elemType = ElemType::Triangle;
        _bftype = BFType::Nedelec;

        geoElem = new RefElem<ElemType::Triangle, BFType::Lagrangian,1>();

        _nG = 3;
        _gW = new double[_nG];
        _gC = new double[_nG * _dim];
        _bf = new double[_nG * _nB * _dim];
        _dbf = new double[_nG * _nB * _dim * _dim];
        _curlbf = new double[_nG * _nB];

        // Quadratures from https://arxiv.org/abs/math/0501496
        //       w_i             x_i                y_i
        _gW[0] = 2.0 / 3.0; _gC[0] = 1.0 / 6.0; _gC[1] = 2.0 / 3.0;
        _gW[1] = 2.0 / 3.0; _gC[2] = 2.0 / 3.0; _gC[3] = 1.0 / 6.0;
        _gW[2] = 2.0 / 3.0; _gC[4] = 1.0 / 6.0; _gC[5] = 1.0 / 6.0;

        double xi, eta;
        for (int iG = 0; iG < _nG; ++iG) {
            xi = _gC[iG * _dim + 0];
            eta = _gC[iG * _dim + 1];

            _bf[iG*_nB*_dim + 0*_dim + 0] = 1.0 - eta; _bf[iG*_nB*_dim + 0*_dim + 1] = xi;
            _bf[iG*_nB*_dim + 1*_dim + 0] = - eta;     _bf[iG*_nB*_dim + 1*_dim + 1] = xi;
            _bf[iG*_nB*_dim + 2*_dim + 0] = - eta;     _bf[iG*_nB*_dim + 2*_dim + 1] = xi - 1.0;

            // dbf/dxi; dbf/deta
            _dbf[iG*_nB*_dim*_dim + 0*_dim*_dim + 0*_dim + 0] = 0.0; _dbf[iG*_nB*_dim*_dim + 0*_dim*_dim + 0*_dim + 1] = -1.0;
            _dbf[iG*_nB*_dim*_dim + 0*_dim*_dim + 1*_dim + 0] = 1.0; _dbf[iG*_nB*_dim*_dim + 0*_dim*_dim + 1*_dim + 1] = 0.0;

            _dbf[iG*_nB*_dim*_dim + 1*_dim*_dim + 0*_dim + 0] = 0.0; _dbf[iG*_nB*_dim*_dim + 1*_dim*_dim + 0*_dim + 1] = -1.0;
            _dbf[iG*_nB*_dim*_dim + 1*_dim*_dim + 1*_dim + 0] = 1.0; _dbf[iG*_nB*_dim*_dim + 1*_dim*_dim + 1*_dim + 1] = 0.0;

            _dbf[iG*_nB*_dim*_dim + 2*_dim*_dim + 0*_dim + 0] = 0.0; _dbf[iG*_nB*_dim*_dim + 2*_dim*_dim + 0*_dim + 1] = -1.0;
            _dbf[iG*_nB*_dim*_dim + 2*_dim*_dim + 1*_dim + 0] = 1.0; _dbf[iG*_nB*_dim*_dim + 2*_dim*_dim + 1*_dim + 1] = 0.0;

            _curlbf[iG*_nB + 0] = 2.0;
            _curlbf[iG*_nB + 1] = 2.0;
            _curlbf[iG*_nB + 2] = 2.0;
        }
    }

    void RefElem<ElemType::Triangle, BFType::Nedelec, 1>::evaluateBFs(const double *refCoords, double *bf) {
        bf[0] = 1.0 - refCoords[1]; bf[1] = refCoords[0];
        bf[2] = - refCoords[1];     bf[3] = refCoords[0];
        bf[4] = - refCoords[1];     bf[5] = refCoords[0] - 1.0;
    }

    void RefElem<ElemType::Triangle, BFType::Nedelec, 1>::evaluateDBFs(const double *refCoords, double *dbf) {
        dbf[0] = 0.0; dbf[1] = -1.0;
        dbf[2] = 1.0; dbf[3] = 0.0;
        dbf[4] = 0.0; dbf[5] = -1.0;
        dbf[6] = 1.0; dbf[7] = 0.0;
        dbf[8] = 0.0; dbf[9] = -1.0;
        dbf[10] = 1.0;dbf[11] = 0.0;
    }

    void RefElem<ElemType::Triangle, BFType::Nedelec, 1>::getJacobian(int iG, const double *nodeCoords, double *J) {
        geoElem->getJacobian(iG, nodeCoords, J);
    }

    void RefElem<ElemType::Triangle, BFType::Nedelec, 1>::getJacobian(const double* dbf, const double* nodeCoords, double* J) {
        geoElem->getJacobian(dbf, nodeCoords, J);
    }

    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


}