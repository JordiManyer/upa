

#include "graph.h"

namespace upa {
    Graph::Graph(int nNodes, int nEdges) {
        _nNodes    = nNodes;
        _nEdges    = nEdges;
        _listSizes = new int[_nNodes+1];
        _adjList   = new int[_nEdges];
    }

    Graph::Graph(int nNodes, int nEdges, const int* sizes) {
        _nNodes    = nNodes;
        _nEdges    = nEdges;
        _listSizes = new int[_nNodes+1];
        _adjList   = new int[_nEdges];

        for (int i = 0; i < _nNodes+1; ++i) _listSizes[i] = sizes[i];
    }

    Graph::Graph(int nNodes, int nEdges, const int* sizes, const int* entries) {
        _nNodes    = nNodes;
        _nEdges    = nEdges;
        _listSizes = new int[_nNodes+1];
        _adjList   = new int[_nEdges];

        for (int i = 0; i < _nNodes+1; ++i) _listSizes[i] = sizes[i];
        for (int i = 0; i < _nEdges; ++i) _adjList[i]   = entries[i];
    }


}

