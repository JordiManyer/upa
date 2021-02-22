
#include <iostream>
#include "structuredMesh.h"
#include "sparse_LIL.h"
#include "sparse_CSR.h"
#include "referenceElement.h"
#include "solver_CG.h"

using namespace std;
using namespace upa;

/***********************************************************************************
 *               TUTORIAL: 2D Laplace with Lagrangian elements
 *
 *     This tutorial solves
 *                   -Lap(u)   = s    ,   (x,y) in [0,1]x[0,1]
 *                   u(x,y)    = 0.0      in  y = 0
 *                   Grad(u)·n = 1.0      in  x = 0, x = 1, y = 1
 *
 *     with analytical solution
 *                   u(x,y) = -s/2 y^2 + sy
 ***********************************************************************************/
int main() {

    double source = 1.0;

    int dim = 2;
    ElemType type = ElemType::Triangle;

    StructuredMesh* mesh = new StructuredMesh();
    mesh->produceCartesian(dim,10,type);

    ReferenceElement* refElem = getReferenceElement(type,BFType::Lagrangian,1);


    int nE = mesh->getNumElements();
    int nN = mesh->getNumNodes();
    int nNbors = mesh->getNumElemNbors();

    int nG = refElem->getNumGaussPoints();
    double* gW = refElem->getGaussWeights();
    double* gC = refElem->getGaussCoords();
    double* bf = refElem->getBasisFunctions(0);
    double* dbf = refElem->getBasisFunctions(1);

    // Count number of nodes in the boundaryDir
    int nBoundary = 0;
    bool boundaryDir[nN]; for (int d = 0; d < nN; ++d) boundaryDir[d] = false;
    bool boundaryNew[nN]; for (int d = 0; d < nN; ++d) boundaryNew[d] = false;
    for (int d = 0; d < nN; ++d) {
        double coords[dim];
        mesh->getNodeCoords(d, coords);

        boundaryDir[d]  = (fabs(coords[1] - 0.0) < 1.e-3);
        boundaryNew[d] = (fabs(coords[0] - 0.0) < 1.e-3) or (fabs(coords[1] - 1.0) < 1.e-3) or (fabs(coords[0] - 1.0) < 1.e-3);
        if (boundaryDir[d]) nBoundary++;
    }

    cout << nBoundary << endl;

    auto sysK = new Sparse_LIL(nN+nBoundary);
    double sysF[nN+nBoundary]; for(int i = 0; i < nN+nBoundary; ++i) sysF[i] = 0.0;
    double Area = 0.0;

    for (int e = 0; e < nE; ++e) { // Loop in elements
        int nodes[nNbors];
        double nodeCoords[nNbors*dim];
        mesh->getElemNodes(e,nodes);
        mesh->getElemCoords(e,nodeCoords);

        double Ke[nNbors*nNbors];
        double fe[nNbors];
        for (int i = 0; i < nNbors; ++i) {
            for (int j = 0; j < nNbors; ++j) Ke[i*nNbors+j] = 0.0;
            fe[i] = 0.0;
        }

        for (int k = 0; k < nG; ++k) { // Loop in Gauss Points

            /// Get data for this gpoint (reference coordinates)
            double wk = gW[k];
            double gCk[dim]; for (int i = 0; i < dim; ++i) gCk[i] = gC[dim*k+i];
            double bfk[nNbors]; for (int i = 0; i < nNbors; ++i) bfk[i] = bf[nNbors*k+i];
            double dbfk[nNbors*dim]; for (int i = 0; i < nNbors*dim; ++i) dbfk[i] = dbf[nNbors*dim*k+i];

            /// Change to physical coordinates
            double J[dim*dim]; // Jacobian
            refElem->getJacobian(k,nodeCoords,J);

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
            Area += dV;

            /// Integration
            // Elemental matrix
            for (int i = 0; i < nNbors; ++i) {
                for (int j = 0; j < nNbors; ++j) {
                    for (int l = 0; l < dim; ++l) Ke[i*nNbors+j] += grad[i*dim+l] * grad[j*dim+l] * dV;
                }
            }
            // Elemental vector
            for (int i = 0; i < nNbors; ++i) {
                fe[i] += bfk[i] * source * dV;
                if (boundaryNew[nodes[i]]) fe[i] += bfk[i] * 0.0 * dV;
            }

        }

        /// Assemble
        for (int i = 0; i < nNbors; ++i) {
            for (int j = 0; j < nNbors; ++j) {
                sysK->assemble(nodes[i], nodes[j], Ke[i*nNbors+j]);
            }
            sysF[nodes[i]] += fe[i];
        }
    }


    /** Apply boundaryDir conditions (using Lagrange multipliers)
     *           | A  K^t |
     * A_final = | K  0   |
     *
     * b_final^t = [b,e]^t
     *
     * with K^t x = e the constraints we want to enforce
     */
    int iB = 0;
    for (int d = 0; d < nN; ++d) {
        if (boundaryDir[d]) {
            sysK->assemble(d,nN+iB,1.0);
            sysK->assemble(nN+iB,d,1.0);
            sysF[nN+iB] = 0.0;
            iB++;
        }
    }


    /// Solve system
    // Convert to CSR
    auto sysK_solve = new Sparse_CSR(sysK);

    cout << "Complete matrix: " << endl;
    for (int i = 0; i < nN+nBoundary; ++i) {
        for (int j = 0; j < nN+nBoundary ; ++j) {
            cout << "    " << sysK_solve->operator()(i,j);
        }
        cout << endl;
    }
    cout << endl;

    cout << "Complete vector: " << endl;
    for (int i = 0; i < nN+nBoundary; ++i) {
        cout << "    " << sysF[i];
    }
    cout << endl;

    // Create solver
    auto CG = new Solver_CG(nN+nBoundary,sysK_solve,sysF);
    CG->setIterations(1000);
    CG->setTolerance(1.e-12);
    CG->setVerbosity(1);

    // Solve
    double x0[nN];
    for (int i = 0; i < nN; ++i) x0[i] = 0.0;
    CG->solve(x0);

    // Output
    double sol[nN+nBoundary];
    CG->getSolution(sol);
    cout << "Converged : " << CG->getConvergence() << endl;
    cout << "Number of iterations: " << CG->getNumIter() << endl;
    cout << "Error2 : " << CG->getError() << endl;

    cout << "x = " << endl;
    for (int i = 0; i < nN; ++i) {
        cout << sol[i] << " , ";
    }
    cout << endl;


    /// Postprocess :: Getting some approximation of the error
    // The analytical solution is u(x) = s/2 ( - x^2 + x)
    double nodalError = 0.0;
    for (int i = 0; i < nN; ++i) {
        double coords[2]; mesh->getNodeCoords(i, coords);
        double analyticalSol = -source/2.0 * coords[1]*coords[1] + source * coords[1];
        nodalError += (sol[i] - analyticalSol)*(sol[i] - analyticalSol);
        cout << analyticalSol << "  ,  ";
    }
    nodalError = sqrt(nodalError);

    cout << endl << "Calculated error of the solution: " << nodalError << endl;

}
