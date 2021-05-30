

#ifndef UPA_STRUCTUREDMESH_H
#define UPA_STRUCTUREDMESH_H

#include <vector>
#include "elementDefinitions.h"

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
        int  getNumEdges() const { return _nEdges;}
        int  getNumNborElems(int e) const { return EEmap[e].size();}
        int  getNumNborNodes(int d) const { return NNmap[d].size();}
        ElemType getElemType() {return _elemType;}

        void getElemNodes(int e, int *dofs);
        void getElemNbors(int e, int *nbors);
        void getElemCoords(int e, double *dofCoords);
        void getElemEdges(int e, int *edges);
        void getNodeCoords(int d, double *dofCoords);
        void getNodeNbors(int d, int *nbors);
        void getEdgeNodes(int e, int *nodes);


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

        // Edge elements
        void produceEdges();

        // Mesh distribution
        void producePartition(int nProcs);

        int _nElems; // Number of elements in the mesh
        int _nNodes; // Number of nodes in the mesh
        int _nNbors; // Number of neighboring nodes per element
        int _dim;    // Physical dimension (1D, 2D, 3D, ...)
        ElemType _elemType;

        int *ENmap;               // Element-Node map, size [nElems,nNbors]
        std::vector<int> *EEmap;  // Element-Element map
        std::vector<int> *NNmap;  // Node-Node map
        double *nodeCoords;       // Node coordinates, size [nNodes,dim]

        // Edge elements
        int _nEdges;
        int *edgeNodes;           // Nodes in each edge [nEdges,2]
        int *edgeMap;             // Element-Edge map [nElems,nNbors]

        // Mesh distribution
        int _nParts;
        int *_epart;
        int *_npart;
    };

}

#endif //UPA_STRUCTUREDMESH_H
