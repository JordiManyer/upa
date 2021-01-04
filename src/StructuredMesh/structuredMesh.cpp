
#include "structuredMesh.h"

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
                    connect[i * _nNbors + 2] = i + nSide + 1 + i/nSide;
                    connect[i * _nNbors + 3] = i + nSide + 2 + i/nSide;
                }
                break;

            default :
                throw std::runtime_error("StructuredMesh: Element Type not recognised.");

        }

    }

}

