
#include "Spline.h"

namespace upa {

    Spline::Spline(int n_in, double *x_in, double *y_in) {
        double h[n_in - 1], alpha[n_in - 1], l[n_in], z[n_in], mu[n_in];

        n = n_in ; x = new double [n];
        a = new double [n]; b = new double [n];
        c = new double [n]; d = new double [n];

        for (int i = 0; i < n; ++i) x[i] = x_in[i];
        for (int i = 0; i < n; ++i) a[i] = y_in[i];
        for (int i = 0; i < n - 1; ++i) h[i] = x[i + 1] - x[i];
        for (int i = 1; i < n - 1; ++i)
            alpha[i] = (3.0 / h[i]) * (a[i + 1] - a[i]) - (3.0 / h[i - 1]) * (a[i] - a[i - 1]);

        l[0] = 1.0;
        mu[0] = 0.0;
        z[0] = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        l[n - 1] = 1.0;
        z[n - 1] = 0.0;
        c[n - 1] = 0.0;

        for (int i = n - 1; i >= 0; --i) {
            c[i] = z[i] - mu[i] * c[i + 1];
            b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
            d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        }

    }


    double Spline::evaluate(double x_eval) {
        int iL, iR, iM;
        double y_eval, h;

        y_eval = 0.0;
        iL = 0;
        iR = n-1;

        // Is x_eval out of the range of our spline?
        if (x_eval < x[iL] or x_eval > x[iR]) throw std::runtime_error("Spline:: Trying to evaluate out of bounds!!");

        // Find interval [x_i , x_{i+1}] containing x_eval
        while (iR != iL) {
            iM = (iL + iR) / 2; // Middle point
            if (x_eval <= x[iM + 1]) iR = iM; // We can shrink the interval by up
            if (x_eval >= x[iM]) iL = iM; // We can shrink the interval by down
        }

        // Evaluate y_eval with the corresponding spline
        h = (x_eval - x[iL]);
        y_eval = a[iL] + b[iL] * h + c[iL] * h * h + d[iL] * h * h * h;
        return y_eval;
    }

}
