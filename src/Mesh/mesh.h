
#ifndef UPA_MESH_H
#define UPA_MESH_H

#include "graph.h"
#include "elementDefinitions.h"
#include "structuredMesh.h"

namespace upa {


/***********************************************************************************************************************
 ***************************                DISTRIBUTED MESH CLASS                           ***************************
 **********************************************************************************************************************/
    class Mesh {

    public:
        Mesh() = default;
        Mesh(StructuredMesh* M);
        ~Mesh() = default;

        /// Connectivity and structure
        int _dim;                 // Physical dimension (1D, 2D, 3D, ...)
        int _nElems;              // Number of elements
        int _nNodes;              // Number of nodes
        int _nNodesI;             // Number of interior nodes
        int _nNodesB;             // Number of boundary nodes
        int _nEdges;              // Number of edges
        int _nNbors;              // Number of neighboring nodes and edges per element
        ElemType _elemType;       // Line, Square, Triangle,

        Graph  *ENmap;           // Element-Node map
        Graph  *EEmap;           // Element-Element map
        Graph  *NNmap;           // Node-Node map
        int    *EEdgeMap;        // Element-Edge map [nElems,nNbors]
        double *nodeCoords;      // Node coordinates, size [nNodes,dim]
        int    *edgeNodes;       // Nodes in each edge [nEdges,2]

        /// MPI and Communication
        int   _rank;              // ID of current proc
        int   _nProcs;            // Number of neighboring procs
        int   *_procs;            // Ranks of neighboring procs
        Graph *_com_nodes;        // Proc-boundNode map
        Graph *_com_edges;        // Proc-boundEdge map


    private:

    };

}


#endif //UPA_MESH_H
