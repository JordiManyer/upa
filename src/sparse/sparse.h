

#ifndef UPA_SPARSE_H
#define UPA_SPARSE_H

#include "sparse_CSR.h"
#include "sparse_LIL.h"
#include "sparse_CSC.h"

// Converts a LIL matrix to CSR
sparse_CSR* LILtoCSR(sparse_LIL* lil);

#endif //UPA_SPARSE_H
