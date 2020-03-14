

#ifndef UPA_SPARSE_CSC_H
#define UPA_SPARSE_CSC_H


class sparse_CSC {

private:
    int n;
    int nonzeros;
    double* rows;
    double* cols;
    double* values;


public:

    sparse_CSC();

    double get(int i, int j);
    void matvec(double* x, double* y);

};


#endif //UPA_SPARSE_CSC_H
