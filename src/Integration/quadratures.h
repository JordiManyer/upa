
#ifndef UPA_QUADRATURES_H
#define UPA_QUADRATURES_H

#include <cmath>
#include <stdexcept>
#include "elementDefinitions.h"

namespace upa {

    void getQuadratureGauss(ElemType type, int degree, int& nG, double* gC, double* gW);

    void getQuadratureGauss_Line(int degree, int& nG, double** gC, double** gW);
    void getQuadratureGauss_Square(int degree, int& nG, double** gC, double** gW);
    void getQuadratureGauss_Triangle(int degree, int& nG, double** gC, double** gW);

}



#endif //UPA_QUADRATURES_H
