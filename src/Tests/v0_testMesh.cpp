

#include <iostream>
#include "structuredMesh.h"

using namespace std;
using namespace upa;


int main() {

    int dim = 2;
    StructuredMesh* mesh = new StructuredMesh();
    mesh->produceCartesian(dim,3,ElemType::Square);

    int nE = mesh->getNumElements();
    int nN = mesh->getNumNodes();
    int nNbors = mesh->getNumElemNbors();

    for (int e = 0; e < nE; ++e) {
        int nodes[nNbors];
        double nodeCoords[nNbors*dim];
        mesh->getElemDOFs(e,nodes);
        mesh->getElemCoords(e,nodeCoords);

        cout << "  -> Element " << e << ":: " << endl;
        for (int i = 0; i < nNbors; ++i) {
            cout << "    -> Node " << nodes[i] << ": (";
            for (int j = 0; j < dim; ++j) cout << nodeCoords[i*dim+j] << ((j != dim-1) ? " , " : "") ;
            cout << ")" << endl;

        }
        cout << endl;

    }

}