

#include "quadratures.h"
#include "myMath.h"


namespace upa {

    void getQuadratureGauss(ElemType type, int degree, int& nG, double** gC, double** gW) {
        switch (type) {
            case ElemType::Line:
                getQuadratureGauss_Line(degree, nG, gC, gW);
                return;

            case ElemType::Square:
                getQuadratureGauss_Square(degree, nG, gC, gW);
                return;

            case ElemType::Triangle:
                getQuadratureGauss_Triangle(degree, nG, gC, gW);
                return;

            default:
                throw std::runtime_error("getQuadratureGauss: Element type not implemented!");
        }
    }



    void getQuadratureGauss_Line(int degree, int& nG, double** gC, double** gW) {
        nG = degree+1;
        *gW = new double[nG];
        *gC = new double[nG];
        GLquad(nG, -1.0, 1.0, *gW, *gC);
    }



    void getQuadratureGauss_Square(int degree, int& nG, double** gC, double** gW) {
        switch (degree) {
            case 1:
                nG = 4;
                *gW = new double[4];
                *gC = new double[8];

                (*gW)[0] = 1.0;           (*gW)[1] = 1.0;          (*gW)[2] = 1.0;          (*gW)[3] = 1.0;           // w
                (*gC)[0] = -1/sqrt(3); (*gC)[2] = 1/sqrt(3);  (*gC)[4] = 1/sqrt(3); (*gC)[6] = -1/sqrt(3); // x
                (*gC)[1] = -1/sqrt(3); (*gC)[3] = -1/sqrt(3); (*gC)[5] = 1/sqrt(3); (*gC)[7] = 1/sqrt(3);  // y
                return;

            case 2:
                nG = 9;
                *gW = new double[9];
                *gC = new double[18];

                (*gW)[0] = 25.0/81.0; (*gW)[1] = 40.0/81.0; (*gW)[2] = 25.0/81.0;
                (*gW)[3] = 40.0/81.0; (*gW)[4] = 64.0/81.0; (*gW)[5] = 40.0/81.0;
                (*gW)[6] = 25.0/81.0; (*gW)[7] = 40.0/81.0; (*gW)[8] = 25.0/81.0;

                (*gC)[0] = - sqrt(3.0/5.0); (*gC)[2] = - sqrt(3.0/5.0); (*gC)[4] = - sqrt(3.0/5.0); // x
                (*gC)[1] = - sqrt(3.0/5.0); (*gC)[3] = 0.0;                (*gC)[5] =   sqrt(3.0/5.0); // y

                (*gC)[6] = 0.0;                (*gC)[8] = 0.0;                (*gC)[10] = 0.0; // x
                (*gC)[7] = - sqrt(3.0/5.0); (*gC)[9] = 0.0;                (*gC)[11] =   sqrt(3.0/5.0); // y

                (*gC)[12] = sqrt(3.0/5.0);   (*gC)[14] = sqrt(3.0/5.0); (*gC)[16] = sqrt(3.0/5.0); // x
                (*gC)[13] = - sqrt(3.0/5.0); (*gC)[15] = 0.0;              (*gC)[17] = sqrt(3.0/5.0); // y
                return;

            default:
                throw std::runtime_error("getQuadratureGauss_square: Quadrature order not implemented!");
        }
    }



    void getQuadratureGauss_Triangle(int degree, int& nG, double** gC, double** gW) {
        switch (degree) {
            case 1:
                nG = 3;
                *gW = new double[3];
                *gC = new double[6];

                //       w_i             x_i                y_i
                (*gW)[0] = 1.0/6.0; (*gC)[0] = 1.0/6.0; (*gC)[1] = 2.0/3.0;
                (*gW)[1] = 1.0/6.0; (*gC)[2] = 2.0/3.0; (*gC)[3] = 1.0/6.0;
                (*gW)[2] = 1.0/6.0; (*gC)[4] = 1.0/6.0; (*gC)[5] = 1.0/6.0;
                return;

            case 2:
                nG = 4;
                *gW = new double[nG];
                *gC = new double[nG*2];

                //       x_i                        y_i                         w_i
                (*gC)[0]  = 0.075031110222608; (*gC)[1]  = 0.280019915499074; (*gW)[0]  = 0.090979309128011;
                (*gC)[2]  = 0.280019915499074; (*gC)[3]  = 0.075031110222608; (*gW)[1]  = 0.090979309128011;
                (*gC)[4]  = 0.178558728263616; (*gC)[5]  = 0.666390246014701; (*gW)[2]  = 0.159020690871988;
                (*gC)[6]  = 0.666390246014701; (*gC)[7]  = 0.178558728263616; (*gW)[3]  = 0.159020690871988;

                return;

            case 3:
                nG = 9;
                *gW = new double[nG];
                *gC = new double[nG*2];

                //       x_i                        y_i                         w_i
                (*gC)[0]  = 0.023931132287081;   (*gC)[1]  = 0.188409405952072; (*gW)[0]  = 0.019396383305959;
                (*gC)[2]  = 0.106170269119576;   (*gC)[3]  = 0.106170269119576; (*gW)[1]  = 0.031034213289535;
                (*gC)[4]  = 0.188409405952072;   (*gC)[5]  = 0.023931132287081; (*gW)[2]  = 0.019396383305959;
                (*gC)[6]  = 0.066554067839164;   (*gC)[7]  = 0.523979067720101; (*gW)[3]  = 0.063678085099885;
                (*gC)[8]  = 0.295266567779633;   (*gC)[9]  = 0.295266567779633; (*gW)[4]  = 0.101884936159816;
                (*gC)[10] = 0.523979067720101;   (*gC)[11] = 0.066554067839165; (*gW)[5]  = 0.063678085099885;
                (*gC)[12] = 0.102717654809626;   (*gC)[13] = 0.808694385677670; (*gW)[6]  = 0.055814420483044;
                (*gC)[14] = 0.455706020243648;   (*gC)[15] = 0.455706020243648; (*gW)[7]  = 0.089303072772871;
                (*gC)[16] = 0.808694385677670;   (*gC)[17] = 0.102717654809626; (*gW)[8]  = 0.055814420483044;
                return;

            default:
                throw std::runtime_error("getQuadratureGauss_triangle: Quadrature order not implemented!");
        }
    }


}