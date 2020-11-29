
#ifndef UPA_INTERPOLATORS_H
#define UPA_INTERPOLATORS_H

const int maxNSplines = 1000;

struct spline {
    int n;
    double a[maxNSplines];
    double b[maxNSplines];
    double c[maxNSplines];
    double d[maxNSplines];
    double x[maxNSplines];
};


void fitCubicNaturalSpline(int n , double* x, double* y, spline* S);
double evalCubicNaturalSpline(double x_eval, spline* S);


#endif //UPA_INTERPOLATORS_H
