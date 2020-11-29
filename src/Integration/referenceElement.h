

#ifndef UPA_REFERENCEELEMENT_H
#define UPA_REFERENCEELEMENT_H

#include <stdexcept>
#include <cmath>

namespace upa {

    enum class ElemType {
        Line, Square
    };

    enum class BFType {
        Lagrangian
    };


    class ReferenceElement {

    public:
        ReferenceElement(ElemType etype, BFType bftype, int bforder);
        ~ReferenceElement() = default;

        // Getters
        int     getNumGaussPoints() {return _nG;}
        double* getGaussWeights() {return _gW;}
        double* getGaussCoords() {return _gC;}
        double* getBasisFunctions(int order) {
            if (order == 0) return _bf;
            if (order == 1) return _dbf;
            else return nullptr;
        }


    private:
        int _dim;
        int _bforder;
        ElemType _elemType;
        BFType _bftype;

        int _nG;      // Number of gauss points
        int _nB;      // Number of basis functions
        double* _gW;  // Gauss weights
        double* _gC;  // Gauss points coordinates
        double* _bf;  // Basis functions evaluated at the gauss points
        double* _dbf; // Derivatives of the basis functions evaluated at the gauss points

        template<ElemType etype, BFType bftype>
        void _fillRefElem(int bforder);

    };



    // Specialisations of _fillRefElem
    template <>
    void ReferenceElement::_fillRefElem<ElemType::Line, BFType::Lagrangian>(int bforder);

    template <>
    void ReferenceElement::_fillRefElem<ElemType::Square, BFType::Lagrangian>(int bforder);

}


#endif //FEMPLAS_REFERENCEELEMENT_H
