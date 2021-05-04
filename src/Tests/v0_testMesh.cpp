

#include <iostream>
#include "structuredMesh.h"

using namespace std;
using namespace upa;


int main() {

    int dim = 2;
    StructuredMesh* mesh = new StructuredMesh();
    mesh->produceCartesian(dim,2,ElemType::Triangle);
    mesh->produceEdges();

    int nE = mesh->getNumElements();
    int nN = mesh->getNumNodes();
    int nNbors = mesh->getNumElemNbors();

    // Print elements
    for (int e = 0; e < nE; ++e) {
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
    }


    // Print nodes
    cout << endl;
    for (int i = 0; i < nN; ++i) {
        int numNbors = mesh->getNumNborNodes(i);
        int nodes[numNbors];
        mesh->getNodeNbors(i,nodes);

        cout << "  -> Node " << i << ":: ";
        for (int j = 0; j < numNbors; ++j) cout << nodes[j] << " , ";
        cout << endl;
    }
    cout << endl;


    // Print edges
    int nEdges = mesh->getNumEdges();
    for (int i = 0; i < nEdges; ++i) {
        int nodes[2]; mesh->getEdgeNodes(i,nodes);
        cout << "  -> Edge " << i << ":: " << nodes[0] << " - " << nodes[1];
        cout << endl;
    }
    cout << endl;

    double p[2] = {0.7,0.9};
    for (int e = 0; e < nE; ++e) if (mesh->isInsideElement(e,p)) cout << "Selected point is inside element " << e << endl;

    int e = mesh->findContainingElem(p);
    cout << "Selected point is inside element " << e << endl;
}