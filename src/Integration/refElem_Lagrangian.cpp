
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
        _dim = 1;
        _nB = 3;
        _bforder = 2;
        _elemType = ElemType::Line;
        _bftype = BFType::Lagrangian;

        _nG = 3;
        _gW = new double[_nG];
        _gC = new double[_nG * _dim];
        _bf = new double[_nG * _nB];
        _dbf = new double[_nG * _nB * _dim];

        _gW[0] = 5.0/9.0;            _gW[1] = 8.0/9.0; _gW[2] = 5.0/9.0;
        _gC[0] = - sqrt(3.0/5.0); _gC[1] = 0.0;     _gC[2] = sqrt(3.0/5.0);

        double xi;
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

    void RefElem<ElemType::Line, BFType::Lagrangian, 2>::evaluateBFs(const double *coords, double *bf) {
        bf[0] = coords[0] * (coords[0] - 1.0) / 2.0;
        bf[1] = 1.0 - coords[0] * coords[0];
        bf[2] = coords[0] * (coords[0] + 1.0) / 2.0;
    }

    void RefElem<ElemType::Line, BFType::Lagrangian, 2>::evaluateDBFs(const double *coords, double *dbf) {
        dbf[0] =  coords[0] - 1.0 / 2.0;
        dbf[1] =  -2.0 * coords[0];
        dbf[2] =  coords[0] - 1.0 / 2.0;
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
        _dim      = 2;
        _nB       = 9;
        _bforder  = 2;
        _elemType = ElemType::Square;
        _bftype   = BFType::Lagrangian;

        _nG  = 9;
        _gW  = new double[_nG];
        _gC  = new double[_nG*_dim];
        _bf  = new double[_nG*_nB];
        _dbf = new double[_nG*_nB*_dim];

        _gW[0] = 25.0/81.0; _gW[1] = 40.0/81.0; _gW[2] = 25.0/81.0;
        _gW[3] = 40.0/81.0; _gW[4] = 64.0/81.0; _gW[5] = 40.0/81.0;
        _gW[6] = 25.0/81.0; _gW[7] = 40.0/81.0; _gW[8] = 25.0/81.0;

        _gC[0] = - sqrt(3.0/5.0); _gC[2] = - sqrt(3.0/5.0); _gC[4] = - sqrt(3.0/5.0); // x
        _gC[1] = - sqrt(3.0/5.0); _gC[3] = 0.0;                _gC[5] =   sqrt(3.0/5.0); // y

        _gC[6] = 0.0;                _gC[8] = 0.0;                _gC[10] = 0.0; // x
        _gC[7] = - sqrt(3.0/5.0); _gC[9] = 0.0;                _gC[11] =   sqrt(3.0/5.0); // y

        _gC[12] = sqrt(3.0/5.0);   _gC[14] = sqrt(3.0/5.0); _gC[16] = sqrt(3.0/5.0); // x
        _gC[13] = - sqrt(3.0/5.0); _gC[15] = 0.0;              _gC[17] = sqrt(3.0/5.0); // y

        double xi,eta;
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

    void RefElem<ElemType::Square, BFType::Lagrangian, 2>::evaluateBFs(const double *coords, double *bf) {
        bf[0] = (coords[0]*(coords[0]-1.0)/2.0) * (coords[1]*(coords[1]-1.0)/2.0);
        bf[1] = (coords[0]*(coords[0]-1.0)/2.0) * (1.0-coords[1]*coords[1]);
        bf[2] = (coords[0]*(coords[0]-1.0)/2.0) * (coords[1]*(coords[1]+1.0)/2.0);

        bf[3] = (1.0-coords[0]*coords[0]) * (coords[1]*(coords[1]-1.0)/2.0);
        bf[4] = (1.0-coords[0]*coords[0]) * (1.0-coords[1]*coords[1]);
        bf[5] = (1.0-coords[0]*coords[0]) * (coords[1]*(coords[1]+1.0)/2.0);

        bf[6] = (coords[0]*(coords[0]+1.0)/2.0) * (coords[1]*(coords[1]-1.0)/2.0);
        bf[7] = (coords[0]*(coords[0]+1.0)/2.0) * (1.0-coords[1]*coords[1]);
        bf[8] = (coords[0]*(coords[0]+1.0)/2.0) * (coords[1]*(coords[1]+1.0)/2.0);
    }

    void RefElem<ElemType::Square, BFType::Lagrangian, 2>::evaluateDBFs(const double *coords, double *dbf) {
        // dN/dxi; dN/deta
        dbf[0] = (coords[0]-1.0/2.0) * (coords[1]*(coords[1]-1.0)/2.0) ;  dbf[1] = (coords[0]*(coords[0]-1.0)/2.0) * (coords[1]-1.0/2.0);
        dbf[2] = (coords[0]-1.0/2.0) * (1.0-coords[1]*coords[1]);         dbf[3] = (coords[0]*(coords[0]-1.0)/2.0) * (- 2.0*coords[1]);
        dbf[4] = (coords[0]-1.0/2.0) * (coords[1]*(coords[1]+1.0)/2.0);   dbf[5] = (coords[0]*(coords[0]-1.0)/2.0) * (coords[1]+1.0/2.0);

        dbf[6] = (-2.0*coords[0]) * (coords[1]*(coords[1]-1.0)/2.0) ;  dbf[7] = (1.0-coords[0]*coords[0]) * (coords[1]-1.0/2.0);
        dbf[8] = (-2.0*coords[0]) * (1.0-coords[1]*coords[1]);         dbf[9] = (1.0-coords[0]*coords[0]) * (- 2.0*coords[1]);
        dbf[10] = (-2.0*coords[0]) * (coords[1]*(coords[1]+1.0)/2.0);   dbf[11] = (1.0-coords[0]*coords[0]) * (coords[1]+1.0/2.0);

        dbf[12] = (coords[0]+1.0/2.0) * (coords[1]*(coords[1]-1.0)/2.0) ;  dbf[13] = (coords[0]*(coords[0]+1.0)/2.0) * (coords[1]-1.0/2.0);
        dbf[14] = (coords[0]+1.0/2.0) * (1.0-coords[1]*coords[1]);         dbf[15] = (coords[0]*(coords[0]+1.0)/2.0) * (- 2.0*coords[1]);
        dbf[16] = (coords[0]+1.0/2.0) * (coords[1]*(coords[1]+1.0)/2.0);   dbf[17] = (coords[0]*(coords[0]+1.0)/2.0) * (coords[1]+1.0/2.0);
    }


    ///**************************************************************************************************************///
    ///**************************************************************************************************************///


}