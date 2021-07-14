
#include "mesh.h"

#include "mpi.h"
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>
#include <iostream>

namespace upa {

    Mesh::Mesh(StructuredMesh *M) {
        int nP;
        MPI_Comm_size(MPI_COMM_WORLD, &nP);
        MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
        std::cout << "Brrrrr" << std::endl;
        if (_rank == 0) std::cerr << nP << " , " << _rank << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);

        if (_rank == 0) M->producePartition(nP);

        int counter;
        int r1, r2, r3, r4, r5; r1 = 10;
        int *rBuffer1, *rBuffer2, *rBuffer3, *rBuffer4, *rBuffer5;
        int EEmap_size, ENmap_size, NNmap_size, com_size;
        if (_rank == 0) {
            int starts[nP], counts[nP], starts_bis[nP];
            int s1, s2, s3, s4, s5;
            int *sBuffer1, *sBuffer2, *sBuffer3, *sBuffer4, *sBuffer5;
            std::set<int> PPmap[nP];
            std::set<int> NPmap[M->_nNodes];

            // Produce Node-Proc and Proc-Proc maps
            for (int iE = 0; iE < M->_nElems; ++iE)
                for (int iN = 0; iN < _nNbors; ++iN) {
                    int node = M->ENmap[iE*_nNbors+iN];
                    if (M->_epart[iE] != M->_npart[node]) M->_npart[node] = -1; // Mark as Boundary node
                    NPmap[node].insert(M->_epart[iE]);
                }

            for (int iN = 0; iN < M->_nNodes; ++iN) if (M->_npart[iN] == -1) {
                for (auto iP = NPmap[iN].begin(); iP != NPmap[iN].end(); ++iP)
                    for (auto iP2 = std::next(iP); iP2 != NPmap[iN].end(); ++iP2) {
                        PPmap[*iP].insert(*iP2);
                        PPmap[*iP2].insert(*iP);
                    }
            }

            // Buffer 1: [_dim, _elemType, _nProcs, _nElems, _nNodesI, _nNodesB, EEmap_size, ENmap_size, NNmap_size, com_size]
            s1 = nP*r1;
            sBuffer1 = new int[s1];
            for (int iP = 0; iP < nP; ++iP) {
                starts[iP] = iP*r1;
                counts[iP] = r1;
                sBuffer1[iP*r1+0] = M->_dim; // _dim
                sBuffer1[iP*r1+1] = static_cast<int>(M->_elemType); // _elemType
                sBuffer1[iP*r1+2] = static_cast<int>(PPmap[iP].size()); // _nProcs
                for (int j = 3; j < r1; ++j) sBuffer1[iP*r1+j] = 0;
            }
            for (int iE = 0; iE < M->_nElems; ++iE) {
                int eProc = M->_epart[iE];
                ++sBuffer1[eProc*r1+3]; // _nElems
                sBuffer1[eProc*r1+6] += static_cast<int>(M->EEmap[iE].size()); // EEmap_size
                sBuffer1[eProc*r1+7] += M->_nNbors; // ENmap_size
                for (int iN = 0; iN < _nNbors; ++iN) {
                    int node = M->ENmap[iE*_nNbors+iN];
                    if (M->_npart[node] == -1) ++sBuffer1[eProc*r1+5]; // _nNodesB
                    else ++sBuffer1[eProc*r1+4]; // _nNodesI
                }
            }
            for (int iN = 0; iN < M->_nNodes; ++iN) {
                for (auto iP : NPmap[iN]) sBuffer1[iP*r1+8] += static_cast<int>(M->NNmap->size()); // NNmap_size
                if (M->_npart[iN] == -1)
                    for (auto iP : NPmap[iN]) sBuffer1[iP*r1+9] += static_cast<int>(NPmap[iN].size())-1; // com_size
            }

            rBuffer1 = new int[r1];
            MPI_Scatterv(sBuffer1, counts, starts, MPI_INT, rBuffer1, r1, MPI_INT, 0, MPI_COMM_WORLD);
            _dim       = rBuffer1[0];
            _elemType  = static_cast<ElemType>(rBuffer1[1]);
            _nProcs    = rBuffer1[2];
            _nElems    = rBuffer1[3];
            _nNodesI   = rBuffer1[4];
            _nNodesB   = rBuffer1[5];
            _nNodes    = rBuffer1[4] + rBuffer1[5];
            EEmap_size = rBuffer1[6];
            ENmap_size = rBuffer1[7];
            NNmap_size = rBuffer1[8];
            com_size   = rBuffer1[9];
            delete[] rBuffer1;

            // Buffer 2: [nborProcs, elems, nodesI, nodesB]
            s2 = 0.0;
            for (int iP = 0; iP < nP; ++iP) {
                starts[iP] = s2;
                counts[iP] = 0.0;
                s2 += sBuffer1[iP*r1+2] + sBuffer1[iP*r1+3] + sBuffer1[iP*r1+4] + sBuffer1[iP*r1+5];
            }
            sBuffer2 = new int[s2];
            for (int iP = 0; iP < nP; ++iP) { // Neighboring partitions
                for (auto iP2 : PPmap[iP]) {
                    sBuffer2[starts[iP] + counts[iP]] = iP2;
                    ++counts[iP];
                }
            }
            for (int iE = 0; iE < M->_nElems; ++iE) { // Elements
                int iP = M->_epart[iE];
                sBuffer2[starts[iP] + counts[iP]] = iE;
                ++counts[iP];
            }
            for (int iN = 0; iN < M->_nNodes; ++iN) if (M->_npart[iN] != -1) { // Interior nodes
                int iP = M->_npart[iN];
                sBuffer2[starts[iP] + counts[iP]] = iN;
                ++counts[iP];
            }
            for (int iN = 0; iN < M->_nNodes; ++iN) if (M->_npart[iN] == -1) { // Boundary nodes
                for (auto iP : NPmap[iN]) {
                    sBuffer2[starts[iP] + counts[iP]] = iN;
                    ++counts[iP];
                }
            }

            r2 = _nProcs + _nElems + _nNodesI + _nNodesB;
            rBuffer2 = new int[r2];
            MPI_Scatterv(sBuffer2, counts, starts, MPI_INT, rBuffer2, r2, MPI_INT, 0, MPI_COMM_WORLD);

            // Buffer 3: [EEmap_sizes, EEmap, ENmap_sizes, ENmap, NNmap_sizes, NNmap]
            s3 = 0;
            for (int iP = 0; iP < nP; ++iP) {
                starts_bis[iP] = starts[iP] + sBuffer1[iP*r1+2];
                starts[iP] = s3;
                counts[iP] = 0.0;
                s3 += 2*sBuffer1[iP*r1+3] + sBuffer1[iP*r1+4] + sBuffer1[iP*r1+5] + sBuffer1[iP*r1+6] + sBuffer1[iP*r1+7] + sBuffer1[iP*r1+8];
            }
            sBuffer3 = new int[s3];
            for (int iP = 0; iP < nP; ++iP) { // For each processor
                for (int iE = starts_bis[iP]; iE < starts_bis[iP]+sBuffer1[iP*r1+3]; ++iE) { // EEmap_sizes
                    sBuffer3[starts[iP] + counts[iP]] = static_cast<int>(M->EEmap[sBuffer2[iE]].size());
                    ++counts[iP];
                }
                for (int iE = starts_bis[iP]; iE < starts_bis[iP]+sBuffer1[iP*r1+3]; ++iE) { // EEmap
                    for (auto iE2 : M->EEmap[sBuffer2[iE]]) {
                        sBuffer3[starts[iP] + counts[iP]] = iE2;
                        ++counts[iP];
                    }
                }
                for (int iE = starts_bis[iP]; iE < starts_bis[iP]+sBuffer1[iP*r1+3]; ++iE) { // ENmap_sizes
                    sBuffer3[starts[iP] + counts[iP]] = M->_nNbors;
                    ++counts[iP];
                }
                for (int iE = starts_bis[iP]; iE < starts_bis[iP]+sBuffer1[iP*r1+3]; ++iE) { // ENmap
                    for (int iN = 0; iN < M->_nNbors; ++iN) {
                        sBuffer3[starts[iP] + counts[iP]] = M->ENmap[sBuffer2[iE]*M->_nNbors+iN];
                        ++counts[iP];
                    }
                }
                starts_bis[iP] = starts_bis[iP] + sBuffer1[iP*r1+3];
                for (int iN = starts_bis[iP]; iN < starts_bis[iP]+sBuffer1[iP*r1+4]+sBuffer1[iP*r1+5]; ++iN) { // NNmap_sizes
                    sBuffer3[starts[iP] + counts[iP]] = static_cast<int>(M->NNmap[sBuffer2[iN]].size());
                    ++counts[iP];
                }
                for (int iN = starts_bis[iP]; iN < starts_bis[iP]+sBuffer1[iP*r1+4]+sBuffer1[iP*r1+5]; ++iN) { // NNmap
                    for (auto iN2 : M->NNmap[sBuffer2[iN]]) {
                        sBuffer3[starts[iP] + counts[iP]] = iN2;
                        ++counts[iP];
                    }
                }
                starts_bis[iP] = starts_bis[iP] + sBuffer1[iP*r1+4];
            }

            r3 =  2*_nElems + _nNodesI + _nNodesB + EEmap_size + ENmap_size + NNmap_size;
            rBuffer3 = new int[r3];
            MPI_Scatterv(sBuffer3, counts, starts, MPI_INT, rBuffer3, r3, MPI_INT, 0, MPI_COMM_WORLD);

            counter = 0;
            EEmap = new Graph(_nElems,EEmap_size);
            ENmap = new Graph(_nElems,ENmap_size);
            NNmap = new Graph(_nNodes,NNmap_size);
            for (int iE = 0; iE < _nElems; ++iE) EEmap->_sizes[iE+1] = EEmap->_sizes[iE] + rBuffer3[counter++];
            for (int iE = 0; iE < EEmap_size; ++iE) EEmap->_nbors[iE] = rBuffer3[counter++];
            for (int iE = 0; iE < _nElems; ++iE) ENmap->_sizes[iE+1] = ENmap->_sizes[iE] + rBuffer3[counter++];
            for (int iN = 0; iN < EEmap_size; ++iN) ENmap->_nbors[iN] = rBuffer3[counter++];
            for (int iN = 0; iN < _nNodes; ++iN) NNmap->_sizes[iN+1] = NNmap->_sizes[iN] + rBuffer3[counter++];
            for (int iN = 0; iN < EEmap_size; ++iN) NNmap->_nbors[iN] = rBuffer3[counter++];
            delete[] rBuffer3;
            delete[] sBuffer3;

            // Buffer 4: [com_size, com_nodes]
            s4 = 0;
            for (int iP = 0; iP < nP; ++iP) {
                starts[iP] = s4;
                counts[iP] = static_cast<int>(PPmap[iP].size());
                s4 += sBuffer1[iP*r1+2] + sBuffer1[iP*r1+9];
            }
            sBuffer4 = new int[s4];
            for (int iP = 0; iP < nP; ++iP) {
                // Count number of boundary nodes shared with each neighboring process
                for (int iN = starts_bis[iP]; iN < starts_bis[iP] + sBuffer1[iP*r1+5]; ++iN) { // Loop through boundary nodes
                    int iN2 = sBuffer2[iN];
                    for (int iP2 : NPmap[iN]) { // iN2 is shared by iP with all neighboring procs iP2
                        if (iP2 != iP) {
                            int it = static_cast<int>(distance(PPmap[iP].begin(), PPmap[iP].find(iP2))); // Get local numbering for iP2
                            ++sBuffer4[starts[iP]+it];
                        }
                    }
                }
                // Get sub-start positions for each neighboring node
                int starts_bisbis[PPmap[iP].size()]; starts_bisbis[0] = starts[iP] + static_cast<int>(PPmap[iP].size());
                for (int iP2 = 0; iP2 < PPmap[iP].size()-1; ++iP2) starts_bisbis[iP2+1] = starts_bisbis[iP2] + sBuffer4[starts[iP]+iP2+1];
                // Get communication lists for iP
                for (int iN = starts_bis[iP]; iN < starts_bis[iP] + sBuffer1[iP*r1+5]; ++iN) { // Loop through boundary nodes
                    int iN2 = sBuffer2[iN];
                    for (int iP2 : NPmap[iN]) { // iN2 is shared by iP with all neighboring procs iP2
                        if (iP2 != iP) {
                            sBuffer4[starts_bisbis[iP]] = iN;
                            ++starts_bisbis[iP];
                        }
                    }
                }
            }

            r4 =  _nProcs + com_size;
            rBuffer4 = new int[r4];
            MPI_Scatterv(sBuffer4, counts, starts, MPI_INT, rBuffer4, r4, MPI_INT, 0, MPI_COMM_WORLD);
            counter = 0;
            _com_nodes = new Graph(_nProcs,com_size);
            for (int iP = 0; iP < _nElems; ++iP) _com_nodes->_sizes[iP+1] = _com_nodes->_sizes[iP] + rBuffer4[counter++];
            for (int iN = 0; iN < com_size; ++iN) _com_nodes->_nbors[iN] = rBuffer4[counter++];

        } else {

            // Buffer 1: [_dim, _elemType, _nProcs, _nElems, _nNodesI, _nNodesB]
            rBuffer1 = new int[r1];
            MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT, rBuffer1, r1, MPI_INT, 0, MPI_COMM_WORLD);
            _dim       = rBuffer1[0];
            _elemType  = static_cast<ElemType>(rBuffer1[1]);
            _nProcs    = rBuffer1[2];
            _nElems    = rBuffer1[3];
            _nNodesI   = rBuffer1[4];
            _nNodesB   = rBuffer1[5];
            _nNodes    = rBuffer1[4] + rBuffer1[5];
            EEmap_size = rBuffer1[6];
            ENmap_size = rBuffer1[7];
            NNmap_size = rBuffer1[8];
            com_size   = rBuffer1[9];
            std::cout << "iP =" << _rank << ", Buffer1= " << rBuffer1[0] << rBuffer1[1] << rBuffer1[2] << rBuffer1[3] << rBuffer1[4] << rBuffer1[5] << std::endl;
            delete[] rBuffer1;

            // Buffer 2: [nborProcs, elems, nodesI, nodesB]
            r2 = _nProcs + _nElems + _nNodesI + _nNodesB;
            rBuffer2 = new int[r2];
            MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT, rBuffer2, r2, MPI_INT, 0, MPI_COMM_WORLD);

            // Buffer 3: [EEmap_sizes, EEmap, ENmap_sizes, ENmap, NNmap_sizes, NNmap]
            r3 =  2*_nElems + _nNodesI + _nNodesB + EEmap_size + ENmap_size + NNmap_size;
            rBuffer3 = new int[r3];
            MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT, rBuffer3, r3, MPI_INT, 0, MPI_COMM_WORLD);

            counter = 0;
            EEmap = new Graph(_nElems,EEmap_size);
            ENmap = new Graph(_nElems,ENmap_size);
            NNmap = new Graph(_nNodes,NNmap_size);
            for (int iE = 0; iE < _nElems; ++iE) EEmap->_sizes[iE+1] = EEmap->_sizes[iE] + rBuffer3[counter++];
            for (int iE = 0; iE < EEmap_size; ++iE) EEmap->_nbors[iE] = rBuffer3[counter++];
            for (int iE = 0; iE < _nElems; ++iE) ENmap->_sizes[iE+1] = ENmap->_sizes[iE] + rBuffer3[counter++];
            for (int iN = 0; iN < EEmap_size; ++iN) ENmap->_nbors[iN] = rBuffer3[counter++];
            for (int iN = 0; iN < _nNodes; ++iN) NNmap->_sizes[iN+1] = NNmap->_sizes[iN] + rBuffer3[counter++];
            for (int iN = 0; iN < EEmap_size; ++iN) NNmap->_nbors[iN] = rBuffer3[counter++];
            delete[] rBuffer3;

            // Buffer 4: [com_size, com_nodes]
            r4 =  _nProcs + com_size;
            rBuffer4 = new int[r4];
            MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT, rBuffer4, r4, MPI_INT, 0, MPI_COMM_WORLD);
            counter = 0;
            _com_nodes = new Graph(_nProcs,com_size);
            for (int iP = 0; iP < _nElems; ++iP) _com_nodes->_sizes[iP+1] = _com_nodes->_sizes[iP] + rBuffer4[counter++];
            for (int iN = 0; iN < com_size; ++iN) _com_nodes->_nbors[iN] = rBuffer4[counter++];

        }



    }

}
