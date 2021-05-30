

#include "graph.h"

namespace upa {
    Graph::Graph(int nNodes, int nEdges) {
        _nN    = nNodes;
        _nE    = nEdges;
        _sizes = new int[_nN+1];
        _nbors   = new int[_nE];
    }

    Graph::Graph(int nNodes, int nEdges, const int* sizes) {
        _nN    = nNodes;
        _nE    = nEdges;
        _sizes = new int[_nN+1];
        _nbors   = new int[_nE];

        for (int i = 0; i < _nN+1; ++i) _sizes[i] = sizes[i];
    }

    Graph::Graph(int nNodes, int nEdges, const int* sizes, const int* entries) {
        _nN    = nNodes;
        _nE    = nEdges;
        _sizes = new int[_nN+1];
        _nbors   = new int[_nE];

        for (int i = 0; i < _nN+1; ++i) _sizes[i] = sizes[i];
        for (int i = 0; i < _nE; ++i) _nbors[i]   = entries[i];
    }


}

