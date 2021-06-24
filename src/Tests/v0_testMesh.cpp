// General libraries
#include <iostream>
#include "mpi.h"

// UPA libraries
#include "structuredMesh.h"
#include "mesh.h"
#include "debugIO.h"

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

    StructuredMesh* mesh;
    if (world_rank == 0) { // Produce and print mesh (serial)

        int dim = 2;
        mesh = new StructuredMesh();
        mesh->produceCartesian(dim,4,ElemType::Square);
        mesh->produceEdges();

        int nE = mesh->getNumElements();
        int nN = mesh->getNumNodes();
        int nNbors = mesh->getNumElemNbors();

        // Print elements
/*        for (int e = 0; e < nE; ++e) {
            int nodes[nNbors];
            double nodeCoords[nNbors*dim];
            mesh->getElemNodes(e,nodes);
            mesh->getElemCoords(e,nodeCoords);

            // Nodes in element
            cout << "  -> Element " << e << ":: " << endl;
            for (int i = 0; i < nNbors; ++i) {
                cout << "    -> Node " << nodes[i] << ": (";
                for (int j = 0; j < dim; ++j) cout << nodeCoords[i*dim+j] << ((j != dim-1) ? " , " : "") ;
                cout << ")" << endl;
            }

            // Neighboring elements
            int nNElems = mesh->getNumNborElems(e);
            int elems[nNElems];
            mesh->getElemNbors(e,elems);
            for (int i = 0; i < nNElems; ++i) {
                cout << "    -> Elem " << elems[i] << endl;
            }
            cout << endl;

            // Edges in element
            int edges[nNbors];
            mesh->getElemEdges(e,edges);
            for (int i = 0; i < nNbors; ++i) {
                cout << "    -> Edge " << edges[i] << endl;
            }
            cout << endl;
        }*/

        // Print nodes
/*        cout << endl;
        for (int i = 0; i < nN; ++i) {
            int numNbors = mesh->getNumNborNodes(i);
            int nodes[numNbors];
            mesh->getNodeNbors(i,nodes);

            cout << "  -> Node " << i << ":: ";
            for (int j = 0; j < numNbors; ++j) cout << nodes[j] << " , ";
            cout << endl;
        }
        cout << endl;*/

        // Print edges
/*        int nEdges = mesh->getNumEdges();
        for (int i = 0; i < nEdges; ++i) {
            int nodes[2]; mesh->getEdgeNodes(i,nodes);
            cout << "  -> Edge " << i << ":: " << nodes[0] << " - " << nodes[1];
            cout << endl;
        }
        cout << endl;*/

        // Partitioning using metis
        cout << "Mesh partitioning: " << endl;
        mesh->producePartition(4);
        printArray(nE,mesh->_epart);
        printArray(nN,mesh->_npart);
        cout << endl;

    }
    if (world_rank == 0) cout << "Distributing the mesh to processors: " << endl;
    //if (world_size > 1) // If run using MPI
    {
        // Sync processes
        MPI_Barrier(MPI_COMM_WORLD);

        // Distribute the mesh between processors
        Mesh* distMesh = new Mesh(mesh);
    }

    MPI_Finalize();
}
