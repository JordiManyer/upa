
#include <iostream>
#include "structuredMesh.h"
#include "sparse_LIL.h"
#include "sparse_CSR.h"
#include "referenceElement.h"
#include "solver_CG.h"
#include "debugIO.h"

using namespace std;
using namespace upa;

/** Basic setup for testing Nédelec FEM elements
 *
 *    We will solve
 *        curl curl E + c E = f     on Omega
 *        E x n = 0                 on delta Omega
 *
 *    which can be used to find diverge-free solutions of the wave equation.
 *    Remark: The selected boundary conditions are imposed as regular Dirichlet BCs by Nedelec elements.
 *    Follows closely https://www.dealii.org/reports/nedelec/nedelec.pdf
 *
 *    To check our implementation, we consider the domain Omega = [0,1]x[0,1]
 *    and the analytical solution
 *                w = x^3 (1 - x)^2 + y^2 (1 - y)^2
 *                E = [dw/dy , -dw/dx] = [4y^3 - 6y^2 + 2y , -5x^4 + 8x^3 - 3x^2],
 *    which fulfills E x n = 0  on delta Omega and is divergence-free by construction.
 *    It can then be calculated that we need to choose as source
 *                f = curl curl E + c E = [ c (4y^3 - 6y^2 +2y) + 12 - 24y , -c (-5x^4 + 8x^3 - 3x^2) + 60x^2 - 48x + 6]
 *
 *    Linear convergence is obtained (as expected).
 */

int main() {

    double source[2]; // f
    double velocity = 1.0; // c

    int dim = 2;
    ElemType type = ElemType::Triangle;

    StructuredMesh* mesh = new StructuredMesh();
    mesh->produceCartesian(dim,8,type);
    mesh->produceEdges();

    ReferenceElement* refElem = getReferenceElement(type,BFType::Nedelec,1);


    int nE   = mesh->getNumElements();
    int nN   = mesh->getNumNodes();
    int nDOF = mesh->getNumEdges();
    int nNbors = mesh->getNumElemNbors();

    int nG = refElem->getNumGaussPoints();
    double* gW = refElem->getGaussWeights();
    double* gC = refElem->getGaussCoords();
    double* bf = refElem->getBasisFunctions(0);
    double* dbf = refElem->getBasisFunctions(1);
    double* curlbf = refElem->getCurlBF();


    // Identify nodes in the Dirichlet boundary.
    int boundaryDir[nN]; for (int d = 0; d < nN; ++d) boundaryDir[d] = 0;
    bool corner[nN]; for (int d = 0; d < nN; ++d) corner[d] = false;
    for (int d = 0; d < nN; ++d) {
        double coords[dim];
        mesh->getNodeCoords(d, coords);

        if (fabs(coords[1] - 0.0) < 1.e-3) boundaryDir[d] = 1;
        if (fabs(coords[0] - 0.0) < 1.e-3) boundaryDir[d] = 2;
        if (fabs(coords[1] - 1.0) < 1.e-3) boundaryDir[d] = 3;
        if (fabs(coords[0] - 1.0) < 1.e-3) boundaryDir[d] = 4;
        corner[d] = ((fabs(coords[1] - 0.0) < 1.e-3) and (fabs(coords[0] - 0.0) < 1.e-3)) or
                    ((fabs(coords[1] - 1.0) < 1.e-3) and (fabs(coords[0] - 0.0) < 1.e-3)) or
                    ((fabs(coords[1] - 1.0) < 1.e-3) and (fabs(coords[0] - 1.0) < 1.e-3)) or
                    ((fabs(coords[1] - 0.0) < 1.e-3) and (fabs(coords[0] - 1.0) < 1.e-3));
    }
    // Identify edges in the Dirichlet boundary
    int nBoundary = 0;
    bool edgesDir[nDOF]; for (int d = 0; d < nN; ++d) edgesDir[d] = false;
    for (int e = 0; e < nDOF; ++e) {
        int nodes[2];
        mesh->getEdgeNodes(e,nodes);
        edgesDir[e] = (boundaryDir[nodes[0]] != 0 and boundaryDir[nodes[1]] != 0 and boundaryDir[nodes[0]] == boundaryDir[nodes[1]]) or
                      (corner[nodes[0]] and boundaryDir[nodes[1]] != 0) or (corner[nodes[1]] and boundaryDir[nodes[0]] != 0);
        if (edgesDir[e]) nBoundary++;
    }

    auto sysK = new Sparse_LIL(nDOF+nBoundary);
    double sysF[nDOF+nBoundary]; for(int i = 0; i < nDOF+nBoundary; ++i) sysF[i] = 0.0;
    double Area = 0.0;

    // Integrate the equations
    for (int e = 0; e < nE; ++e) { // Loop in elements
        int nodes[nNbors], edges[nNbors], edgeSigns[nNbors];
        double nodeCoords[nNbors*dim];
        mesh->getElemNodes(e,nodes);
        mesh->getElemCoords(e,nodeCoords);
        mesh->getElemEdges(e,edges);

        for (int i = 0; i < nNbors; ++i) { // Specific for triangles. Careful!
            if ( i != nNbors-1 and nodes[i] < nodes[i+1]) edgeSigns[i] = 1;
            else if ( i == nNbors-1 and nodes[i] < nodes[0]) edgeSigns[i] = 1;
            else edgeSigns[i] = -1;
        }

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
            double bfk[nNbors*dim]; for (int i = 0; i < nNbors*dim; ++i) bfk[i] = bf[nNbors*dim*k+i];
            double dbfk[nNbors*dim*dim]; for (int i = 0; i < nNbors*dim*dim; ++i) dbfk[i] = dbf[nNbors*dim*dim*k+i];
            double curlbfk[nNbors]; for (int i = 0; i < nNbors; ++i) curlbfk[i] = curlbf[nNbors*k+i];

            /// Change to physical coordinates
            double J[dim*dim]; // Jacobian, dx/dxi
            refElem->getJacobian(k,nodeCoords,J);

            double Jinv[dim*dim]; inverse(dim, J, Jinv);
            double detJ = det(dim,J);

            double bfk_xy[nNbors*dim]; // Nedelec basis functions in physical coordinates.
            double curlbfk_xy[nNbors]; // Curl in physical coordinates.
            for (int l = 0; l < nNbors; ++l) {
                curlbfk_xy[l] = curlbfk[l]/detJ;
                for (int i = 0; i < dim; ++i) {
                    bfk_xy[l*dim + i] = 0.0;
                    for (int j = 0; j < dim; ++j) {
                        bfk_xy[l*dim + i] += Jinv[i*dim + j] * bfk[l*dim + j];
                    }
                }
            }

            double dV = wk * detJ;
            Area += dV;

            // Calculate source f = curl curl E + c E
            double pCoords[dim]; // Physical coordinates (x,y) of the Gauss point in this element
            refElem->getPhysicalCoords(k,nodeCoords,pCoords);
            double x = pCoords[0]; double y = pCoords[1];
            source[0] =  velocity*(4*y*y*y - 6*y*y+ 2*y) + (12 - 24*y);
            source[1] = -velocity*(-5*x*x*x*x+ 8*x*x*x - 3*x*x) + (60*x*x - 48*x + 6);

            /// Integration
            // Elemental matrix
            for (int i = 0; i < nNbors; ++i) {
                for (int j = 0; j < nNbors; ++j) {
                    Ke[i*nNbors+j] += curlbfk_xy[i] * curlbfk_xy[j] * dV;
                    for (int l = 0; l < dim; ++l) Ke[i*nNbors+j] += velocity * bfk_xy[i*dim+l] * bfk_xy[j*dim+l] * dV;
                }
            }
            // Elemental vector
            for (int i = 0; i < nNbors; ++i) {
                for (int l = 0; l < dim; ++l) fe[i] += bfk_xy[i*dim+l] * source[l] * dV;
            }
        }

        /// Assemble
        for (int i = 0; i < nNbors; ++i) {
            for (int j = 0; j < nNbors; ++j) {
                sysK->assemble(edges[i], edges[j], edgeSigns[i] * edgeSigns[j] * Ke[i*nNbors+j]);
            }
            sysF[edges[i]] += edgeSigns[i] * fe[i];
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
    for (int e = 0; e < nDOF; ++e) {
        if (edgesDir[e]) {
            sysK->assemble(e,nDOF+iB,1.0);
            sysK->assemble(nDOF+iB,e,1.0);
            sysF[nDOF+iB] = 0.0;
            iB++;
        }
    }

    /// Solve system
    // Convert to CSR
    auto sysK_solve = new Sparse_CSR(sysK);

    cout << "Complete matrix: " << endl;
    for (int i = 0; i < nDOF+nBoundary; ++i) {
        for (int j = 0; j < nDOF+nBoundary ; ++j) {
            cout << "    " << sysK_solve->operator()(i,j);
        }
        cout << endl;
    }
    cout << endl;

    cout << "Complete vector: " << endl;
    for (int i = 0; i < nDOF+nBoundary; ++i) {
        cout << "    " << sysF[i];
    }
    cout << endl;

    // Create solver
    auto CG = new Solver_CG(nDOF+nBoundary,sysK_solve,sysF);
    CG->setIterations(1000);
    CG->setTolerance(1.e-12);
    CG->setVerbosity(1);

    // Solve
    double x0[nDOF+nBoundary];
    for (int i = 0; i < nDOF+nBoundary; ++i) x0[i] = 0.0;
    CG->solve(x0);

    // Output
    double sol[nDOF+nBoundary];
    CG->getSolution(sol);
    cout << "Converged : " << CG->getConvergence() << endl;
    cout << "Number of iterations: " << CG->getNumIter() << endl;
    cout << "Error2 : " << CG->getError() << endl;

    cout << "x = " << endl;
    for (int i = 0; i < nDOF+nBoundary; ++i) {
        cout << sol[i] << " , ";
    }
    cout << endl;



    /// L2 Error of the solution, i.e Err = integral { |sol - analyticalSolution|^2 dOmega }
    cout << endl;
    cout << "Comparing with analytical solution: " << endl;

    double Error = 0.0;
    Area = 0.0;
    for (int e = 0; e < nE; ++e) { // Loop in elements
        int nodes[nNbors], edges[nNbors], edgeSigns[nNbors];
        double nodeCoords[nNbors*dim];
        mesh->getElemNodes(e,nodes);
        mesh->getElemCoords(e,nodeCoords);
        mesh->getElemEdges(e,edges);

        for (int i = 0; i < nNbors; ++i) { // Specific for triangles. Careful!
            if ( i != nNbors-1 and nodes[i] < nodes[i+1]) edgeSigns[i] = 1;
            else if ( i == nNbors-1 and nodes[i] < nodes[0]) edgeSigns[i] = 1;
            else edgeSigns[i] = -1;
        }

        for (int k = 0; k < nG; ++k) { // Loop in Gauss Points
            double dofk[nNbors];
            for (int i = 0; i < nNbors; ++i) dofk[i] = sol[edges[i]];

            double wk = gW[k];
            double bfk[nNbors*dim]; for (int i = 0; i < nNbors*dim; ++i) bfk[i] = bf[nNbors*dim*k+i];
            double J[dim*dim]; refElem->getJacobian(k,nodeCoords,J);
            double Jinv[dim*dim]; inverse(dim,J,Jinv);
            double detJ = det(dim,J);

            /// TODO: This should be done as part of interpolateSolution() in Nedelec!
            double approxE[2], auxE[2];
            for (int i = 0; i < dim; ++i) {
                auxE[i] = 0.0;
                for (int j = 0; j < nNbors; ++j) {
                    auxE[i] += edgeSigns[j] * bfk[j*dim+i] * dofk[j];
                }
            }
            for (int i = 0; i < dim; ++i) { // Needs to be multiplied by the Jacobian as part of Piola Transformation.
                approxE[i] = 0.0;
                for (int j = 0; j < dim; ++j) {
                    approxE[i] += Jinv[i*dim + j] * auxE[j];
                }
            }

            // Get exact solution
            double pCoords[dim];
            refElem->getPhysicalCoords(k,nodeCoords,pCoords);
            double analE[2];
            double x = pCoords[0]; double y = pCoords[1];
            analE[0] =  (4*y*y*y - 6*y*y+ 2*y);
            analE[1] = (-5*x*x*x*x+ 8*x*x*x - 3*x*x);

            double dV = wk * detJ;
            Error += ((analE[0] - approxE[0]) * (analE[0] - approxE[0]) + (analE[1] - approxE[1]) * (analE[1] - approxE[1])) * dV;
            Area += dV;
        }
    }

    Error = sqrt(Error);
    cout << "L2 error :: " << Error << endl;
    cout << endl;

}