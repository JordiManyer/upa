#ifndef UPA_GRAPH_H
#define UPA_GRAPH_H

namespace upa {

    /**
     ** @class Directed Graph, represented by a storage-optimized adjacency list.
     **/
    class Graph {
    public:

        // Constructors and destructors
        Graph(int nNodes, int nEdges);
        Graph(int nNodes, int nEdges, const int* sizes);
        Graph(int nNodes, int nEdges, const int* sizes, const int* entries);
        ~Graph() = default;

        // Read / Write operators
        int operator()(int i, int j) const { return _nbors[_sizes[i] + j]; }
        int &operator()(int i, int j) { return _nbors[_sizes[i] + j]; }

        int _nN;     // Number of nodes in the graph
        int _nE;     // Number of directed edges in the graph
        int *_sizes; // Indexes where each sublist starts (size _nN+1)
        int *_nbors; // Concatenated lists of neighbors for each node (size _nE)
    };

}

#endif //UPA_GRAPH_H
