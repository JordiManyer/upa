
#include <iostream>
#include "structuredMesh.h"
#include "sparse_LIL.h"
#include "sparse_CSR.h"
#include "referenceElement.h"
#include "solver_CG.h"

using namespace std;
using namespace upa;


int main() {

    double source = 1.0;

    int dim = 1;
    ElemType type = ElemType::Line;

    auto* mesh = new StructuredMesh();
    mesh->produceCartesian(dim,12,type);

    auto* refElem = new ReferenceElement(type,BFType::Lagrangian,1);


    int nE = mesh->getNumElements();
    int nN = mesh->getNumNodes();
    int nNbors = mesh->getNumElemNbors();

    int nG = refElem->getNumGaussPoints();
    double* gW = refElem->getGaussWeights();
    double* gC = refElem->getGaussCoords();
    double* bf = refElem->getBasisFunctions(0);
    double* dbf = refElem->getBasisFunctions(1);

    auto sysK = new Sparse_LIL(nN);
    double sysF[nN]; for(int i = 0; i < nN; ++i) sysF[i] = 0.0;

    for (int e = 0; e < nE; ++e) { // Loop in elements
        int nodes[nNbors];
        double nodeCoords[nNbors*dim];
        mesh->getElemNodes(e,nodes);
        mesh->getElemCoords(e,nodeCoords);


        for (int k = 0; k < nG; ++k) { // Loop in Gauss Points

            /// Get data for this gpoint (reference coordinates)
            double wk = gW[k];
            double gCk[dim]; for (int i = 0; i < dim; ++i) gCk[i] = gC[dim*k+i];
            double bfk[nNbors]; for (int i = 0; i < nNbors; ++i) bfk[i] = bf[nNbors*k+i];
            double dbfk[nNbors*dim]; for (int i = 0; i < nNbors*dim; ++i) dbfk[i] = dbf[nNbors*dim*k+i];

            /// Change to physical coordinates
            double J[dim*dim]; // Jacobian
            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j) {
                    J[i * dim + j] = 0.0;
                    for (int l = 0; l < nNbors; ++l)
                        J[i * dim + j] += dbfk[dim * l + i] * nodeCoords[dim * l + j];
                }

            double Jinv[dim*dim]; inverse(dim, J, Jinv);

            double grad[nNbors*dim];
            for (int l = 0; l < nNbors; ++l) {
                for (int i = 0; i < dim; ++i) {
                    grad[l*dim + i] = 0.0;
                    for (int j = 0; j < dim; ++j) {
                        grad[l*dim + i] += Jinv[i*dim + j] * dbfk[l*dim + j];
                    }
                }
            }

            double dV = wk * det(dim, J);

            /// Integration
            double Ke[nNbors*nNbors];
            double fe[nNbors];

            // Elemental matrix
            for (int i = 0; i < nNbors; ++i) {
                for (int j = 0; j < nNbors; ++j) {
                    Ke[i*nNbors+j] = 0.0;
                    for (int l = 0; l < dim; ++l) Ke[i*nNbors+j] += grad[i*dim+l] * grad[j*dim+l];
                    Ke[i*nNbors+j] *= dV;
                }
            }
            // Elemental vector
            for (int i = 0; i < nNbors; ++i) {
                fe[i]  = bfk[i] * source * dV;
            }

            /// Assemble
            for (int i = 0; i < nNbors; ++i) {
                for (int j = 0; j < nNbors; ++j) {
                    sysK->assemble(nodes[i], nodes[j], Ke[i*nNbors+j]);
                }
                sysF[nodes[i]] += fe[i];
            }
        }
    }

    /// Dirichlet BCs (Elimination)
    bool select[nN]; int nSelect = 0;
    for (int d = 0; d < nN; ++d) {
        double coords[dim];
        mesh->getNodeCoords(d, coords);

        bool boundary = (fabs(coords[0] - 0.0) < 1.e-3) or (fabs(coords[0] - 1.0) < 1.e-3);
        select[d] = not boundary;
        if (not boundary) nSelect++;
    }

    // Apply BCs on the vector
    double sysF_solve[nSelect];
    int k = 0; double value = 0.0;
    for (int i = 0; i < nN; ++i) {
        if (select[i]) {
            sysF_solve[k] = sysF[i];
            k++;
        }
    }

    cout << "Selected rows: " << endl;
    for (int i = 0; i < nN; ++i) cout << select[i] << " ";
    cout << endl;

    // Convert matrix to CSR + impose elimination BCs
    auto sysK_solve = new Sparse_CSR(sysK,select);

    cout << "Complete matrix: " << endl;
    for (int i = 0; i < sysK_solve->n; ++i) {
        for (int j = 0; j < sysK_solve->n ; ++j) {
            cout << "    " << sysK_solve->operator()(i,j);
        }
        cout << endl;
    }
    cout << endl;

    cout << "Complete Vector: " << endl;
    for (int i = 0; i < nSelect; ++i) cout << "    " << sysF_solve[i];
    cout << endl;


    /// Solve system
    // Create solver
    auto CG = new Solver_CG(nSelect,sysK_solve,sysF_solve);
    CG->setTolerance(0.000001);
    CG->setVerbosity(1);

    // Solve
    double x0[nSelect];
    for (int i = 0; i < nN; ++i) x0[i] = 1.0;
    CG->solve(x0);

    // Output
    double sol[nSelect];
    CG->getSolution(sol);
    cout << "Converged : " << CG->getConvergence() << endl;
    cout << "Number of iterations: " << CG->getNumIter() << endl;
    cout << "Error2 : " << CG->getError() << endl;

    cout << "x = " << endl;
    k = 0;
    for (int i = 0; i < nN; ++i) {
        if (select[i]) {
            cout << sol[k] << " , "; k++;
        } else cout << value << " , ";
    }
    cout << endl;

}

