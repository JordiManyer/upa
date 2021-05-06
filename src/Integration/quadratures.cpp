

#include "quadratures.h"


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
        switch (degree) {
            case 1:
                nG = 2;
                *gW = new double[2];
                *gC = new double[2];

                (*gW)[0] = 1.0; (*gW)[1] = 1.0;
                (*gC)[0] = -1 / sqrt(3); (*gC)[1] = 1 / sqrt(3);
                return;

            case 2:
                nG = 3;
                *gW = new double[3];
                *gC = new double[3];

                (*gW)[0] = 5.0/9.0;            (*gW)[1] = 8.0/9.0; (*gW)[2] = 5.0/9.0;
                (*gC)[0] = - sqrt(3.0/5.0); (*gC)[1] = 0.0;     (*gC)[2] = sqrt(3.0/5.0);
                return;

            default:
                throw std::runtime_error("getQuadratureGauss: Quadrature order not implemented!");
        }
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
                throw std::runtime_error("getQuadratureGauss: Quadrature order not implemented!");
        }
    }



    void getQuadratureGauss_Triangle(int degree, int& nG, double** gC, double** gW) {
        switch (degree) {
            case 1:
                nG = 3;
                *gW = new double[3];
                *gC = new double[6];

                // Quadratures from https://arxiv.org/abs/math/0501496
                //       w_i             x_i                y_i
                (*gW)[0] = 1.0/6.0; (*gC)[0] = 1.0/6.0; (*gC)[1] = 2.0/3.0;
                (*gW)[1] = 1.0/6.0; (*gC)[2] = 2.0/3.0; (*gC)[3] = 1.0/6.0;
                (*gW)[2] = 1.0/6.0; (*gC)[4] = 1.0/6.0; (*gC)[5] = 1.0/6.0;
                return;

            case 2:
                nG = 9;
                *gW = new double[9];
                *gC = new double[18];

                // Quadratures from https://arxiv.org/abs/math/0501496
                //       w_i                        x_i                         y_i
                (*gW)[0] = 0.2199034873106; (*gC)[0]  = 0.0915762135098; (*gC)[1]  = 0.0915762135098;
                (*gW)[1] = 0.2199034873106; (*gC)[2]  = 0.8168475729805; (*gC)[3]  = 0.0915762135098;
                (*gW)[2] = 0.2199034873106; (*gC)[4]  = 0.0915762135098; (*gC)[5]  = 0.8168475729805;
                (*gW)[3] = 0.4467631793560; (*gC)[6]  = 0.1081030181681; (*gC)[7]  = 0.4459484909160;
                (*gW)[4] = 0.4467631793560; (*gC)[8]  = 0.4459484909160; (*gC)[9]  = 0.1081030181681;
                (*gW)[5] = 0.4467631793560; (*gC)[10] = 0.4459484909160; (*gC)[11] = 0.4459484909160;
                return;

            default:
                throw std::runtime_error("getQuadratureGauss: Quadrature order not implemented!");
        }
    }


}