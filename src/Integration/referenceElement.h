

#ifndef UPA_REFERENCEELEMENT_H
#define UPA_REFERENCEELEMENT_H

#include <stdexcept>
#include <cmath>

namespace upa {

    /* Element types. Parenthesis -> Node IDs; Brackets -> Node coordinates
     *
     *   -Line:
     *          (0) ------ (1)
     *          [-1]        [1]
     *
     *   -Square (equivalent to Line x Line):
     *        [-1,1]      [1,1]
     *          (3) ------ (2)
     *           |          |
     *           |          |
     *          (0) ------ (1)
     *        [-1,-1]     [1,-1]
     *
     */
    enum class ElemType {
        Line, Square
    };

    // Basis function types
    enum class BFType {
        Lagrangian // Polynomial basis such that Ni(xj) = delta_ij
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

        /** @brief Fills the ReferenceElement object with all needed information
         *
         * @tparam etype  - Element type
         * @tparam bftype - Basis function type
         * @param bforder - Element order
         */
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
