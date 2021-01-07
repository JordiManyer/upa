
#include "structuredMesh.h"
#include <iostream>

namespace upa {

    void StructuredMesh::getElemDOFs(int e, int *dofs) {
        if (e > _nElems) throw std::runtime_error("StructuredMesh: Element ID out of range.");
        else for (int i = 0; i < _nNbors; ++i) dofs[i] = connect[e*_nNbors + i];
    }

    void StructuredMesh::getElemCoords(int e, double *dofCoords) {
        if (e > _nElems) throw std::runtime_error("StructuredMesh: Element ID out of range.");
        else {
            for (int i = 0; i < _nNbors; ++i) {
                int nodeID = connect[e * _nNbors + i];
                for (int j = 0; j < _dim; ++j) dofCoords[i*_dim + j] = nodeCoords[nodeID*_dim + j];
            }
        }
    }

    void StructuredMesh::getNodeCoords(int d, double *dofCoords) {
        if (d > _nNodes) throw std::runtime_error("StructuredMesh: Node ID out of range.");
        else for (int i = 0; i < _dim; ++i) dofCoords[i] = nodeCoords[d*_dim + i];
    }

    void StructuredMesh::produceCartesian(int dim, int nSide, ElemType type) {
        _dim = dim;
        _nElems = ipow(nSide,dim);
        _nNodes = ipow(nSide+1,dim);
        _elemType = type;

        // Calculate node coordinates
        nodeCoords = new double[_nNodes*_dim];
        int iSide[_dim]; for (int i = 0; i < _dim; ++i) iSide[i] = 0;
        for (int i = 0; i < _nNodes; ++i) {
            for (int j = 0; j < _dim; ++j)
                nodeCoords[i * _dim + j] = double(iSide[j]) / double(nSide);

            ++iSide[0];
            for (int j = 0; j < dim-1; ++j)
                if(iSide[j] == nSide+1) {
                    iSide[j] = 0; ++iSide[j+1];
                }
        }

        // Create connectivity matrix
        switch (_elemType) {
            case ElemType::Line:
                if (dim != 1) throw std::runtime_error("StructuredMesh: ElemType::Line needs dim = 1");
                _nNbors = 2;
                connect = new int[_nElems*_nNbors];

                for (int i = 0; i < _nElems; ++i) {
                    connect[i*_nNbors + 0] = i;
                    connect[i*_nNbors + 1] = i+1;
                }
                break;

            case ElemType::Square:
                if (dim != 2) throw std::runtime_error("StructuredMesh: ElemType::Square needs dim = 2");
                _nNbors = 4;
                connect = new int[_nElems*_nNbors];

                for (int i = 0; i < _nElems; ++i) {
                    connect[i * _nNbors + 0] = i + i/nSide;
                    connect[i * _nNbors + 1] = i + 1 + i/nSide;
                    connect[i * _nNbors + 2] = i + nSide + 2 + i/nSide;
                    connect[i * _nNbors + 3] = i + nSide + 1 + i/nSide;
                }
                break;

            default :
                throw std::runtime_error("StructuredMesh: Element Type not recognised.");

        }

    }

    bool StructuredMesh::isInsideElement(int elem, double *coords) {
        if (_dim == 1) {
            // A point is inside a line if the sign of the product of the signed distances to the vertices is negative,
            // i.e. if the point has one vertex at each side.
            return (nodeCoords[connect[elem*_nNbors+0]] - coords[0]) * (nodeCoords[connect[elem*_nNbors+1]] - coords[0]) < 0;
        } else if (_dim == 2) {
            // A point P is inside a convex polygon, given by vertices A_0,A_1, ... A_n,
            // if and only if all of the a_i = x_{i+1} y_i - x_i y_{i+1} have the same sign
            // where [x_i , y_i] = A_i - P
            double x1[_dim], x2[_dim], a1, a2;
            int nodes[_nNbors]; getElemDOFs(elem,nodes);

            for (int i = 0; i < _dim; ++i) x2[i] = nodeCoords[nodes[_nNbors-1]*_dim+i] - coords[i];
            for (int j = 0; j < _nNbors; ++j) {
                for (int i = 0; i < _dim; ++i) x1[i] = nodeCoords[nodes[j]*_dim+i] - coords[i];
                a1 = x2[0] * x1[1] - x1[0] * x2[1];
/*                std::cout << elem << " " << a1 << std::endl;
                std::cout << nodeCoords[nodes[j]*_dim+0] << "," << nodeCoords[nodes[j]*_dim+1] << std::endl;
                std::cout << x1[0] << "," << x1[1] << std::endl;
                std::cout << x2[0] << "," << x2[1] << std::endl << std::endl;*/
                if (j != 0 and a1 * a2 < 0.0) return false;
                for (int i = 0; i < _dim; ++i) x2[i] = x1[i]; a2 = a1;
            }
            return true;
        } else if (_dim == 3) {
            /// TODO : Implement interior check for 3D elements
            throw std::runtime_error("StructuredMesh: isInsideElement still not implemented for 3D.");
        }

        return false;
    }

    void StructuredMesh::getElemBarycenter(int elem, double *coords) {
        int nodes[_nNbors]; getElemDOFs(elem,nodes);

        for (int i = 0; i < _dim; ++i) coords[i] = 0.0;
        for (int j = 0; j < _nNbors; ++j)
            for (int i = 0; i < _dim; ++i) coords[i] += nodeCoords[nodes[j]*_dim+i];
        for (int i = 0; i < _dim; ++i) coords[i] /= _nNbors;
    }

}

