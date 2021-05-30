
#include <iostream>
#include "mpi.h"
#include "graph.h"

using namespace std;
using namespace upa;

int main() {
    MPI_Init(nullptr, nullptr);
    cout << endl;

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //cout << endl << world_rank << endl;


    if (world_rank == 0) {
        int nNodes = 3;
        int nEdges = 6;
        int sizes[4] = {0,2,4,6};
        int adjList[6] = {1,2,0,2,0,1};

        Graph & G = *(new Graph(nNodes,nEdges,sizes,adjList));
        cout << G(0,0) << endl;
    }


    MPI_Finalize();
}

