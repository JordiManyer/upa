
#include <iostream>
#include <cmath>
#include "structuredMesh.h"
#include "sparse_LIL.h"
#include "sparse_CSR.h"
#include "referenceElement.h"
#include "solver_CG.h"
#include "debugIO.h"
#include "myMath.h"

using namespace std;
using namespace upa;

/***********************************************************************************
 *               TUTORIAL: 1D Laplace with Lagrangian elements
 *
 *     This tutorial solves
 *                   -d^2u/dx^2 = s,    x in [0,1]
 *                   u(0) = u(1) = 0.0
 ***********************************************************************************/
int main() {

    /// Parameters
    int dim = 1;                         // Dimension of the problem
    ElemType etype = ElemType::Line;     // Type of element
    BFType bftype  = BFType::Lagrangian; // Type of basis function
    int bforder    = 1;                  // Order of basis function
    int nElems = 5;                     // Number of elements in the mesh

    double source = 1.0; // Source ('s' in the equation), here constant in the domain


    /// Creation of the mesh and the reference element
    auto* mesh = new StructuredMesh();
    mesh->produceCartesian(dim,nElems,etype);
    auto* refElem = getReferenceElement(etype,bftype,bforder);

    /// Get information from the mesh and the reference element (not necessary, here for clarity)
    int nE = mesh->getNumElements();                    // Number of elements in the mesh
    int nN = mesh->getNumNodes();                       // Number of nodes in the mesh
    int nNbors = mesh->getNumElemNbors();               // Number of nodes per element

    int nG = refElem->getNumGaussPoints();              // Number of gauss points per element
    double* gW = refElem->getGaussWeights();            // Gauss weights, size (nG)
    double* gC = refElem->getGaussCoords();             // Gauss points coordinates, size (nG,dim)
    double* bf = refElem->getBasisFunctions(0);   // Basis functions evaluated at gauss points, size (nG)
    double* dbf = refElem->getBasisFunctions(1);  // Basis function derivatives evaluated at gauss points, size (nG,dim)

    /// Creation of the system matrix and vector, which are going to be filled in the integration loop
    auto sysK = new Sparse_LIL(nN);
    double sysF[nN]; for(int i = 0; i < nN; ++i) sysF[i] = 0.0;

    /** Integration loop
     *   We want to do
     *
     *      Kij = integral {  dNi/dx * dNj/dx  dx}
     *      Fi  = integral {  source * Ni dx }
     *
     *      We will integrate this in each element, using Gauss quadrature. This results in a two loop algorithm:
     *           - The first loop goes through all elements in the mesh and does the inetgral in that element.
     *             Using the quadratures, the integral results in
     *               sum over k { wk * (dNi/dx)(xk) * (dNj/dx)(xk) }
     *             where wk,xk are the gauss weights and gauss points
     *           - Therefore we need a second loop, which does the sum mentioned above in each element.
     *
     */
    double area = 0.0;
    for (int e = 0; e < nE; ++e) { // Loop in elements

        // Create elemental matrix and vector.
        // These will be filled and then assembled into the final matrix and vector.
        double Ke[nNbors*nNbors]; // Elemental matrix
        double fe[nNbors];        // Elemental vector
        for (int i = 0; i < nNbors; ++i) {
            fe[i] = 0.0;
            for (int j = 0; j < nNbors; ++j) Ke[i * nNbors + j] = 0.0;
        }

        // Get the nodes IDs and node coordinates for this element
        int nodes[nNbors];
        double nodeCoords[nNbors*dim];
        mesh->getElemNodes(e,nodes);
        mesh->getElemCoords(e,nodeCoords);

        /// Integration inside this element
        for (int k = 0; k < nG; ++k) { // Loop in Gauss Points

            /// Get data for this gauss point (reference coordinates)
            double wk = gW[k];
            double gCk[dim]; for (int i = 0; i < dim; ++i) gCk[i] = gC[dim*k+i];
            double bfk[nNbors]; for (int i = 0; i < nNbors; ++i) bfk[i] = bf[nNbors*k+i];
            double dbfk[nNbors*dim]; for (int i = 0; i < nNbors*dim; ++i) dbfk[i] = dbf[nNbors*dim*k+i];

            /// Change from reference coordinates to physical coordinates
            double J[dim*dim]; // Jacobian of the change, J_lm = d x_l / d Eta_m
            refElem->getJacobian(k,nodeCoords,J);

            double Jinv[dim*dim]; inverse(dim, J, Jinv); // Invert the Jacobian

            // The values of the basis functions do not change (bfk), but the derivatives do.
            // So we need to compute the new derivatives dN/dx = dEta/dx  * dN/dEta
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
            area += dV;

            /// Integration
            for (int i = 0; i < nNbors; ++i) {
                for (int j = 0; j < nNbors; ++j) {
                    for (int l = 0; l < dim; ++l) Ke[i*nNbors+j] += grad[i*dim+l] * grad[j*dim+l] * dV;
                }
            }
            for (int i = 0; i < nNbors; ++i) {
                fe[i]  += bfk[i] * source * dV;
            }

        } // End of gauss point loop

        /// Assemble into global system
        for (int i = 0; i < nNbors; ++i) {
            for (int j = 0; j < nNbors; ++j) {
                sysK->assemble(nodes[i], nodes[j], Ke[i*nNbors+j]);
            }
            sysF[nodes[i]] += fe[i];
        }

    } // End of element loop

    /** After the loop, the system matrix and vector are completely assembled. However, we need to impose
     ** Dirichlet BCs. In this example we will use elimination. We simply eliminate the corresponding rows and
     ** columns of the system matrix and vector, by imposing the known values.
     */

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

    // Convert matrix to CSR + impose elimination BCs
    auto sysK_solve = new Sparse_CSR(sysK,select);


    /// Print final matrix and vector
    cout << "Complete matrix: " << endl;
    for (int i = 0; i < sysK_solve->n; ++i) {
        for (int j = 0; j < sysK_solve->n ; ++j) {
            cout << "    " << sysK_solve->operator()(i,j);
        }
        cout << endl;
    }
    cout << endl;

    cout << "Complete Vector: " << endl;
    printArray(nSelect,sysF_solve);
    cout << endl;


    /// Solve system
    // Create solver
    auto CG = new Solver_CG(nSelect,sysK_solve,sysF_solve);
    CG->setTolerance(0.000001);
    CG->setVerbosity(1);

    // Solve
    double x0[nSelect];
    for (int i = 0; i < nSelect; ++i) x0[i] = 1.0;
    CG->solve(x0);

    // Output
    double sol[nN];
    double aux[nSelect];
    CG->getSolution(aux);
    k = 0;
    for (int i = 0; i < nN; ++i) {
        if (select[i]) {
            sol[i] = aux[k]; k++;
        } else sol[i] = value;
    }

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
        double x; mesh->getNodeCoords(i, &x);
        double analyticalSol = source/2.0 * ( -x*x + x);
        nodalError += (sol[i] - analyticalSol)*(sol[i] - analyticalSol);
        cout << analyticalSol << "  ,  ";
    }
    nodalError = sqrt(nodalError);

    cout << endl << "Calculated error of the solution: " << nodalError << endl;

}

