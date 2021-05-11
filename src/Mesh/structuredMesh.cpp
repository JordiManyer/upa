
#include "structuredMesh.h"

namespace upa {

    void StructuredMesh::getElemNodes(int e, int *dofs) {
        if (e > _nElems) throw std::runtime_error("Mesh: Element ID out of range.");
        else for (int i = 0; i < _nNbors; ++i) dofs[i] = ENmap[e*_nNbors + i];
    }

    void StructuredMesh::getElemCoords(int e, double *dofCoords) {
        if (e > _nElems) throw std::runtime_error("Mesh: Element ID out of range.");
        else {
            for (int i = 0; i < _nNbors; ++i) {
                int nodeID = ENmap[e * _nNbors + i];
                for (int j = 0; j < _dim; ++j) dofCoords[i*_dim + j] = nodeCoords[nodeID*_dim + j];
            }
        }
    }

    void StructuredMesh::getNodeCoords(int d, double *dofCoords) {
        if (d > _nNodes) throw std::runtime_error("Mesh: Node ID out of range.");
        else for (int i = 0; i < _dim; ++i) dofCoords[i] = nodeCoords[d*_dim + i];
    }

    void StructuredMesh::getElemNbors(int e, int *nbors) {
        if (e > _nElems) throw std::runtime_error("Mesh: Elem ID out of range.");
        else for (int i = 0; i < EEmap[e].size(); ++i) nbors[i] = EEmap[e][i];
    }

    void StructuredMesh::getNodeNbors(int d, int *nbors) {
        if (d > _nNodes) throw std::runtime_error("Mesh: Node ID out of range.");
        else for (int i = 0; i < NNmap[d].size(); ++i) nbors[i] = NNmap[d][i];
    }

    void StructuredMesh::getElemEdges(int e, int *edges) {
        if (e > _nElems) throw std::runtime_error("Mesh: Elem ID out of range.");
        else for (int i = 0; i < _nNbors; ++i) edges[i] = edgeMap[e*_nNbors+i];
    }

    void StructuredMesh::getEdgeNodes(int e, int *nodes) {
        if (e > _nEdges) throw std::runtime_error("Mesh: Edge ID out of range.");
        else for (int i = 0; i < 2; ++i) nodes[i] = edgeNodes[e*2+i];
    }

    void StructuredMesh::produceCartesian(int dim, int nSide, ElemType type) {
        _dim = dim;
        _elemType = type;

        // Calculate node coordinates
        _nNodes = ipow(nSide+1,dim);
        nodeCoords = new double[_nNodes*_dim];
        int iSide[_dim]; for (int i = 0; i < _dim; ++i) iSide[i] = 0;
        for (int i = 0; i < _nNodes; ++i) {
            for (int j = 0; j < _dim; ++j)
                nodeCoords[i * _dim + j] = double(iSide[j]) / double(nSide);

            ++iSide[0];
            for (int j = 0; j < dim-1; ++j)
                if(iSide[j] == nSide+1) {
                    iSide[j] = 0; ++iSide[j+1];
                }
        }

        // Create connectivity matrix
        switch (_elemType) {
            case ElemType::Line:
                if (dim != 1) throw std::runtime_error("Mesh: ElemType::Line needs dim = 1");
                _nElems = ipow(nSide,dim);
                _nNbors = 2;
                ENmap = new int[_nElems*_nNbors];
                EEmap = new std::vector<int>[_nElems];
                NNmap = new std::vector<int>[_nNodes];

                for (int i = 0; i < _nElems; ++i) {
                    ENmap[i*_nNbors + 0] = i;
                    ENmap[i*_nNbors + 1] = i+1;
                }
                for (int i = 0; i < _nElems; ++i) {
                    if (i != 0) EEmap[i].push_back(i-1);
                    if (i != _nElems-1) EEmap[i].push_back(i+1);
                }
                for (int i = 0; i < _nNodes; ++i) {
                    if (i != 0) NNmap[i].push_back(i-1);
                    if (i != _nNodes-1) NNmap[i].push_back(i+1);
                }
                break;

            case ElemType::Square:
                if (dim != 2) throw std::runtime_error("Mesh: ElemType::Square needs dim = 2");
                _nElems = ipow(nSide,dim);
                _nNbors = 4;
                ENmap = new int[_nElems*_nNbors];
                EEmap = new std::vector<int>[_nElems];
                NNmap = new std::vector<int>[_nNodes];

                for (int i = 0; i < _nElems; ++i) {
                    ENmap[i * _nNbors + 0] = i + i/nSide;
                    ENmap[i * _nNbors + 1] = i + 1 + i/nSide;
                    ENmap[i * _nNbors + 2] = i + nSide + 2 + i/nSide;
                    ENmap[i * _nNbors + 3] = i + nSide + 1 + i/nSide;
                }
                for (int i = 0; i < nSide; ++i) {
                    for (int j = 0; j < nSide; ++j) {
                        if (i != 0)       EEmap[i*nSide+j].push_back((i-1)*nSide + j);
                        if (i != nSide-1) EEmap[i*nSide+j].push_back((i+1)*nSide + j);
                        if (j != 0)       EEmap[i*nSide+j].push_back(i*nSide + j-1);
                        if (j != nSide-1) EEmap[i*nSide+j].push_back(i*nSide + j+1);
                    }
                }
                for (int i = 0; i < nSide+1; ++i) {
                    for (int j = 0; j < nSide+1; ++j) {
                        if (i != 0)     NNmap[i*(nSide+1)+j].push_back((i-1)*(nSide+1) + j);
                        if (j != 0)     NNmap[i*(nSide+1)+j].push_back(i*(nSide+1) + j-1);
                        if (j != nSide) NNmap[i*(nSide+1)+j].push_back(i*(nSide+1) + j+1);
                        if (i != nSide) NNmap[i*(nSide+1)+j].push_back((i+1)*(nSide+1) + j);
                    }
                }
                break;

            case ElemType::Triangle:
                if (dim != 2) throw std::runtime_error("Mesh: ElemType::Triangle needs dim = 2");
                _nElems = 2*ipow(nSide,dim);
                _nNbors = 3;
                ENmap = new int[_nElems*_nNbors];
                EEmap = new std::vector<int>[_nElems];
                NNmap = new std::vector<int>[_nNodes];

                for (int i = 0; i < _nElems/2; ++i) { // For each square, we do bottom triangle then top triangle
                    ENmap[(2*i) * _nNbors + 0] = i + i/nSide;
                    ENmap[(2*i) * _nNbors + 1] = i + 1 + i/nSide;
                    ENmap[(2*i) * _nNbors + 2] = i + nSide + 1 + i/nSide;

                    ENmap[(2*i+1) * _nNbors + 0] = i + nSide + 2 + i/nSide;
                    ENmap[(2*i+1) * _nNbors + 1] = i + nSide + 1 + i/nSide;
                    ENmap[(2*i+1) * _nNbors + 2] = i + 1 + i/nSide;
                }
                for (int i = 0; i < nSide; ++i) { // For each square, we do bottom triangle then top triangle
                    for (int j = 0; j < nSide; ++j) {
                        if (i != 0)       EEmap[i*2*nSide+2*j].push_back((i-1)*2*nSide+2*j+1);
                        if (j != 0)       EEmap[i*2*nSide+2*j].push_back(i*2*nSide+2*j-1);
                        EEmap[i*2*nSide+2*j].push_back(i*2*nSide+2*j+1);

                        EEmap[i*2*nSide+2*j+1].push_back(i*2*nSide+2*j);
                        if (j != nSide-1) EEmap[i*2*nSide+2*j+1].push_back(i*2*nSide+2*j+2);
                        if (i != nSide-1) EEmap[i*2*nSide+2*j+1].push_back((i+1)*2*nSide+2*j);
                    }
                }
                for (int i = 0; i < nSide+1; ++i) {
                    for (int j = 0; j < nSide+1; ++j) {
                        if (i != 0)                NNmap[i*(nSide+1)+j].push_back((i-1)*(nSide+1) + j);
                        if (i != 0 and j != nSide) NNmap[i*(nSide+1)+j].push_back((i-1)*(nSide+1) + j+1);
                        if (j != 0)                NNmap[i*(nSide+1)+j].push_back(i*(nSide+1) + j-1);
                        if (j != nSide)            NNmap[i*(nSide+1)+j].push_back(i*(nSide+1) + j+1);
                        if (i != nSide and j != 0) NNmap[i*(nSide+1)+j].push_back((i+1)*(nSide+1) + j-1);
                        if (i != nSide)            NNmap[i*(nSide+1)+j].push_back((i+1)*(nSide+1) + j);
                    }
                }
                break;

            default :
                throw std::runtime_error("Mesh: Element Type not recognised.");

        }

    }

    bool StructuredMesh::isInsideElement(int elem, double *coords) {
        if (_dim == 1) {
            // A point is inside a line if the sign of the product of the signed distances to the vertices is negative,
            // i.e. if the point has one vertex at each side.
            return (nodeCoords[ENmap[elem*_nNbors+0]] - coords[0]) * (nodeCoords[ENmap[elem*_nNbors+1]] - coords[0]) < 0;
        } else if (_dim == 2) {
            // A point P is inside a convex polygon, given by vertices A_0,A_1, ... A_n,
            // if and only if all of the a_i = x_{i+1} y_i - x_i y_{i+1} have the same sign
            // where [x_i , y_i] = A_i - P
            double x1[_dim], x2[_dim], a1, a2;
            int nodes[_nNbors]; getElemNodes(elem,nodes);

            for (int i = 0; i < _dim; ++i) x2[i] = nodeCoords[nodes[_nNbors-1]*_dim+i] - coords[i];
            for (int j = 0; j < _nNbors; ++j) {
                for (int i = 0; i < _dim; ++i) x1[i] = nodeCoords[nodes[j]*_dim+i] - coords[i];
                a1 = x2[0] * x1[1] - x1[0] * x2[1];

                if (j != 0 and a1 * a2 < 0.0) return false;
                for (int i = 0; i < _dim; ++i) x2[i] = x1[i]; a2 = a1;
            }
            return true;
        } else if (_dim == 3) {
            /// TODO : Implement interior check for 3D elements
            throw std::runtime_error("Mesh: isInsideElement still not implemented for 3D.");
        }

        return false;
    }

    void StructuredMesh::getElemBarycenter(int elem, double *coords) {
        int nodes[_nNbors]; getElemNodes(elem,nodes);

        for (int i = 0; i < _dim; ++i) coords[i] = 0.0;
        for (int j = 0; j < _nNbors; ++j)
            for (int i = 0; i < _dim; ++i) coords[i] += nodeCoords[nodes[j]*_dim+i];
        for (int i = 0; i < _dim; ++i) coords[i] /= _nNbors;
    }

    int StructuredMesh::findContainingElem(double *coords, int e0) {
        double dist; double elemCoords[_dim];
        std::pair<double,int> p;
        std::priority_queue<std::pair<double,int>> eq;

        // Push initial search elements
        if (e0 == -1) e0 = _nElems/2;
        dist = 0.0;
        getElemBarycenter(e0,elemCoords);
        for (int i = 0; i < _dim; ++i) dist += fabs(coords[i]-elemCoords[i]);
        eq.push(std::make_pair(-dist,e0));

        // Dijkstra-like search algorithm
        while (not eq.empty()) {
            p = eq.top(); eq.pop();
            if (isInsideElement(p.second,coords)) return p.second; // If found, return element
            for (int e : EEmap[p.second]) {                        // Else, put neighbor elems in queue
                dist = 0.0;
                getElemBarycenter(e,elemCoords);
                for (int i = 0; i < _dim; ++i) dist += fabs(coords[i]-elemCoords[i]);
                eq.push(std::make_pair(-dist,e));
            }
        }
        return -1;
    }

    void StructuredMesh::produceEdges() {
        if (_dim < 2) throw std::runtime_error("Mesh: produceEdges() only for dim > 1.");

        // Count edges
        _nEdges = 0;
        for (int i = 0; i < _nNodes; ++i) _nEdges += NNmap[i].size();
        _nEdges /= 2;

        // Find edge ordering
        edgeNodes = new int[_nEdges*2];
        int k = 0; int indexes[_nNodes+1]; indexes[0] = 0;
        for (int d = 0; d < _nNodes; ++d) {
            for (auto i : NNmap[d])
                if (i > d) {
                    edgeNodes[k*2+0] = d;
                    edgeNodes[k*2+1] = i;
                    k++;
                }
            indexes[d+1] = k;
        }

        // Create Element-Edge map
        edgeMap = new int[_nElems*_nNbors];
        for (int e = 0; e < _nElems; ++e) {
            // Get list of edges for this element
            int edges[_nNbors*2];
            if (_elemType == ElemType::Square) {
                edges[0] = ENmap[e*_nNbors+0]; edges[1] = ENmap[e*_nNbors+1];
                edges[2] = ENmap[e*_nNbors+1]; edges[3] = ENmap[e*_nNbors+2];
                edges[4] = ENmap[e*_nNbors+2]; edges[5] = ENmap[e*_nNbors+3];
                edges[6] = ENmap[e*_nNbors+3]; edges[7] = ENmap[e*_nNbors+0];
            } else if (_elemType == ElemType::Triangle) {
                edges[0] = ENmap[e*_nNbors+0]; edges[1] = ENmap[e*_nNbors+1];
                edges[2] = ENmap[e*_nNbors+1]; edges[3] = ENmap[e*_nNbors+2];
                edges[4] = ENmap[e*_nNbors+2]; edges[5] = ENmap[e*_nNbors+0];
            }

            // For each element edge, look for its numbering.
            for (int i = 0; i < _nNbors; ++i) {
                if (edges[i*2+0] > edges[i*2+1]) { // We want lower-indexed node first.
                    int swap = edges[i*2+0]; edges[i*2+0] = edges[i*2+1]; edges[i*2+1] = swap;
                }
                k = indexes[edges[i*2+0]];
                while (k < indexes[edges[i*2+0]+1] and edgeNodes[k*2+1] != edges[i*2+1]) ++k;
                edgeMap[e*_nNbors+i] = k;
            }
        }
    }

}

