

#ifndef UPA_SPARSE_CSR_H
#define UPA_SPARSE_CSR_H


class sparse_CSR {

private:
    int n;
    int nonzeros;
    double* rows;
    double* cols;
    double* values;


public:

    sparse_CSR();

    double get(int i, int j);
    void matvec(double* x, double* y);

};


#endif //UPA_SPARSE_CSR_H
