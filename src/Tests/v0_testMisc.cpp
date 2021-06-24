
#include <iostream>
#include "mpi.h"
#include "graph.h"
#include "structuredMesh.h"
#include "mesh.h"

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

    auto SMesh = new StructuredMesh();
    SMesh->produceCartesian(2,2,ElemType::Square);

    auto DMesh = new Mesh(SMesh);

    MPI_Finalize();
}

