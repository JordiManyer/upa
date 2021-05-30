
#ifndef UPA_MESH_H
#define UPA_MESH_H

#include "mpi.h"
#include "graph.h"
#include "elementDefinitions.h"

namespace upa {


/***********************************************************************************************************************
 ***************************                DISTRIBUTED MESH CLASS                           ***************************
 **********************************************************************************************************************/
    class Mesh {

    public:
        Mesh() = default;
        ~Mesh() = default;

        /// Connectivity and structure
        int _dim;                 // Physical dimension (1D, 2D, 3D, ...)
        int _nElems;              // Number of elements
        int _nNodes;              // Number of nodes
        int _nEdges;              // Number of edges
        int _nNbors;              // Number of neighboring nodes and edges per element
        ElemType _elemType;       // Line, Square, Triangle,

        int    *_ENmap;           // Element-Node map, size [nElems,nNbors]
        Graph  *_EEmap;           // Element-Element map
        Graph  *_NNmap;           // Node-Node map
        int    *_EEdgeMap;        // Element-Edge map [nElems,nNbors]
        double *_nodeCoords;      // Node coordinates, size [nNodes,dim]
        int    *_edgeNodes;       // Nodes in each edge [nEdges,2]

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
