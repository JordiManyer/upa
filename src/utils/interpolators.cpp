
#include "interpolators.h"
#include <math.h>
#include <iostream>


// Given a set of points (x_i, y_i) i = 0, ... n
// returns a set of n-1 cubic natural splines defined in intervals [x_i, x_{i+1}]
// A natural spline is defined in [x_i, x_{i+1}] and given by the formula
//            S_i(x) = a_i + b_i (x - x_i) + c_i (x - x_i)^2 + d_i (x - x_i)^3
// for i = 0, ... n-1
void fitCubicNaturalSpline(int n , double* x, double* y, spline* S) {

    // Pointers to Spline coefficients
    double *a , *b , *c , *d;
    a = S->a; b = S->b; c = S->c; d = S->d;
    // Auxiliar vars
    double h[n-1] , alpha[n-1], l[n], z[n], mu[n];

    S->n = n-1;
    for (int i = 0; i < n ; ++i)   S->x[i] = x[i];
    for (int i = 0; i < n ; ++i)   a[i] = y[i];
    for (int i = 0; i < n-1 ; ++i) h[i] = x[i+1] - x[i];
    for (int i = 1; i < n-1 ; ++i) alpha[i] = (3.0/h[i]) * (a[i+1]-a[i]) - (3.0/h[i-1]) * (a[i]-a[i-1]);

    l[0] = 1.0; mu[0] = 0.0; z[0] = 0.0;
    for (int i = 1; i < n-1 ; ++i) {
        l[i] = 2.0 * (x[i+1] - x[i-1]) - h[i-1]*mu[i-1];
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i];
    }
    l[n-1] = 1.0; z[n-1] = 0.0; c[n-1] = 0.0;

    for (int i = n-1 ; i >= 0 ; --i) {
        c[i] = z[i] - mu[i]*c[i+1];
        b[i] = (a[i+1] - a[i]) / h[i] - h[i] * (c[i+1] + 2.0*c[i]) / 3.0;
        d[i] = (c[i+1] - c[i]) / (3.0 * h[i]);
    }

}

// Given a point x_eval and a set of splines S, returns the fitted value y_eval = S(x_eval)
double evalCubicNaturalSpline(double x_eval, spline* S) {
    int iL, iR, iM;
    double y_eval, h;

    double* x = S->x;

    y_eval = 0.0;
    iL = 0; iR = S->n;

    // Is x_eval out of the range of our spline?
    if (x_eval < x[iL] or x_eval > x[iR]) std::cout << "Error: Out of bounds in point " << x_eval << std::endl;

    // Find interval [x_i , x_{i+1}] containing x_eval
    while ( iR != iL ) {
        iM = (iL + iR)/2.0; // Middle point
        if (x_eval <= x[iM+1]) iR = iM; // We can shrink the interval by up
        if (x_eval >= x[iM]) iL = iM; // We can shrink the interval by down
    }

    // Evaluate y_eval with the corresponding spline
    h = (x_eval - x[iL]);
    y_eval = S->a[iL] + \
             S->b[iL] * h + \
             S->c[iL] * h * h + \
             S->d[iL] * h * h * h ;
    return y_eval;
}

