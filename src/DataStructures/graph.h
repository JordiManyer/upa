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
        int operator()(int i, int j) const { return _adjList[_listSizes[i] + j]; }
        int &operator()(int i, int j) { return _adjList[_listSizes[i] + j]; }

        int _nNodes;
        int _nEdges;
        int *_listSizes;
        int *_adjList;
    };

}

#endif //UPA_GRAPH_H
