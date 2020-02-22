

#ifndef UPA_FITTERS_H
#define UPA_FITTERS_H


struct data1D {
    int n;
    double *x;
    double *y;
};

struct data2D {
    int n;
    double *x1;
    double *x2;
    double *y;
};


void PLS1D(data1D* data_in, int degree, double* c);
void MLS1D(data1D* data_in, int degree, double rho, data1D* data_out);


#endif //UPA_FITTERS_H
