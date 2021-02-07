

#ifndef UPA_STRUCTUREDMESH_H
#define UPA_STRUCTUREDMESH_H

#include <stdexcept>
#include <vector>
#include <queue>
#include <utility>
#include <iostream>

#include "myMath.h"
#include "referenceElement.h"

namespace upa {

    class StructuredMesh {

    public:
        StructuredMesh() = default;
        ~StructuredMesh() = default;

        // Mesh creation
        void produceCartesian(int dim, int nSide, ElemType type);

        // Getters
        int  getNumElements() const { return _nElems; }
        int  getNumNodes() const { return _nNodes;}
        int  getNumElemNbors() const { return _nNbors;}
        int  getNumNborElems(int e) const { return EEmap[e].size();}
        ElemType getElemType() {return _elemType;}

        void getElemNodes(int e, int *dofs);
        void getElemNbors(int e, int *nbors);
        void getElemCoords(int e, double *dofCoords);
        void getNodeCoords(int d, double *dofCoords);

        // Auxiliar functions

        /** @brief Returns true if the point given by 'coords' is inside the element 'elem'
         */
        bool isInsideElement(int elem, double* coords);

        /** @brief Fills 'coords' with the baricenter coordinates of the element 'elem'
         */
        void getElemBarycenter(int elem, double* coords);

        /** @brief Returns the element ID containing the point 'coords'
         */
        int findContainingElem(double* coords, int e0 = -1);

    private:
        int _nElems; // Number of elements in the mesh
        int _nNodes; // Number of nodes in the mesh
        int _nNbors; // Number of neighboring nodes per element
        int _dim;    // Physical dimension (1D, 2D, 3D, ...)
        ElemType _elemType;

        int *ENmap;               // Element-Node map, size [nElems,nNbors]
        std::vector<int> *EEmap;  // Element-Element map
        double *nodeCoords;       // Node coordinates, size [nNodes,dim]

    };

}

#endif //UPA_STRUCTUREDMESH_H
