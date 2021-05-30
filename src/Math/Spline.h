
#ifndef UPA_SPLINE_H
#define UPA_SPLINE_H

namespace upa {

    /** @brief Cubic natural spline
    *   Given a set of points (x_i, y_i) i = 0, ... n
    *   returns a set of n-1 cubic natural splines defined in intervals [x_i, x_{i+1}]
    *   A natural spline is defined in [x_i, x_{i+1}] and given by the formula
    *              S_i(x) = a_i + b_i (x - x_i) + c_i (x - x_i)^2 + d_i (x - x_i)^3
    *   for i = 0, ... n-1
    */
    class Spline {

    public:
        Spline(int n_in, const double* x_in, const double* y_in);
        ~Spline() = default;

        /** @brief Evaluate spline
         *
         * @param x_eval Evaluation point
         * @return y_eval = Spline(x_eval)
         */
        double evaluate(double x_eval);

    private:
        int n;
        double* a;
        double* b;
        double* c;
        double* d;
        double* x;
    };

}

#endif //UPA_SPLINE_H
