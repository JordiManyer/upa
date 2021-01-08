

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
        mesh->getElemNodes(e,nodes);
        mesh->getElemCoords(e,nodeCoords);

        cout << "  -> Element " << e << ":: " << endl;
        for (int i = 0; i < nNbors; ++i) {
            cout << "    -> Node " << nodes[i] << ": (";
            for (int j = 0; j < dim; ++j) cout << nodeCoords[i*dim+j] << ((j != dim-1) ? " , " : "") ;
            cout << ")" << endl;
        }

        int nNElems = mesh->getNumNborElems(e);
        int elems[e];
        mesh->getElemNbors(e,elems);
        for (int i = 0; i < nNElems; ++i) {
            cout << "    -> Elem " << elems[i] << endl;
        }
        cout << endl;
    }

    double p[2] = {0.7,0.8};
    for (int e = 0; e < nE; ++e) if (mesh->isInsideElement(e,p)) cout << "Selected point is inside element " << e << endl;

    int e = mesh->findContainingElem(p);
    cout << "Selected point is inside element " << e << endl;
}