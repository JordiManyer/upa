#include "utils/fitters.h"
#include "utils/interpolators.h"
#include <math.h>
#include <iostream>


void testFitters();
void testInterpolators();
double testFunc(double x);
double testDFunc(double x);

int main() {
    testFitters();
    testInterpolators();
    return 0;
}


void testFitters() {
    double a, b, sum;
    a = 0.0; b = 10.0;

    // Create data to be fitted
    int n = 100;
    double x[n], y[n], yfit[n];
    for (int i = 0; i < n; ++i) {
        x[i] = a + i * (b - a) / (n - 1);
        y[i] = testFunc(x[i]);
    }
    data1D data_in, data_out;
    data_in.n = n; data_in.x = x; data_in.y = y;
    data_out.n = n; data_out.x = x; data_out.y = yfit;

    // Polynomial least squares
    std::cout << "Testing Polynomial least squares fitting: " << std::endl;
    int degree = 2; double c[degree+1];
    PLS1D(&data_in, degree, c);
    std::cout << "Degree of fitting = " << degree << std::endl;
    std::cout << "PLS coeffs c0, c1, c2 : " << c[0] << " " << c[1] << " " << c[2] << std::endl;

    double ePLS = 0.0;
    for (int i = 0; i < n; ++i) {
        sum = 0.0;
        for (int j = 0; j < degree+1; ++j) sum += c[j] * pow(data_in.x[i],j);
        ePLS += sqrt(fabs(data_in.y[i]*data_in.y[i] - sum*sum));
    }
    std::cout << "Error of fitting: Total =  " << ePLS  << " Per point = " << ePLS/n << std::endl << std::endl;

    // Moving least squares
    std::cout << "Testing moving least squares fitting: " << std::endl;
    std::cout << "Degree of fitting = " << degree << std::endl;

    double rho = 3.2;
    MLS1D( &data_in, degree, rho, &data_out);

    double eMLS = 0.0;
    for (int i = 0; i < n; ++i) {
        eMLS += sqrt(fabs(data_in.y[i]*data_in.y[i] - data_out.y[i]*data_out.y[i]));
    }
    std::cout << "Error of fitting: Total =  " << eMLS  << " Per point = " << eMLS/n << std::endl << std::endl;
}


void testInterpolators () {
    double a, b;
    a = 0.0; b = 10.0;

    // Create data to be fitted
    int n = 100;
    double x[n], y[n];
    for (int i = 0; i < n; ++i) {
        x[i] = a + i * (b - a) / (n - 1);
        y[i] = testFunc(x[i]);
    }

    // Test cubic natural splines
    std::cout << "Testing Cubic Natural Splines: " << std::endl;
    // Create spline
    spline S;
    fitCubicNaturalSpline(n , x, y, &S);
    // Evaluate spline
    double y_eval[n];
    for (int i = 0; i < n; ++i) y_eval[i] = evalCubicNaturalSpline(x[i], &S);
    // Evaluate spline error
    double eSPL = 0.0;
    for (int i = 0; i < n; ++i) eSPL += sqrt(fabs( y[i]*y[i] - y_eval[i]*y_eval[i] ));
    std::cout << "Error of interpolation: Total =  " << eSPL  << " Per point = " << eSPL/n << std::endl << std::endl;

}



double testFunc(double x) {
    return 1.0 + 2.0*x - 1.0*x*x;
}

double testDFunc(double x) {
    return 2.0 - 2.0*x;
}