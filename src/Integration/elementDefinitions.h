

#ifndef UPA_ELEMENTDEFINITIONS_H
#define UPA_ELEMENTDEFINITIONS_H

namespace upa {

    /***************    ELEMENT TYPES    ********************
     * Parenthesis -> Node IDs; Brackets -> Node coordinates
     *
     *   - Line:
     *          (0) ------ (1)
     *          [-1]        [1]
     *
     *   - Square (equivalent to Line x Line):
     *        [-1,1]      [1,1]
     *          (3) ------ (2)
     *           |          |
     *           |          |
     *          (0) ------ (1)
     *        [-1,-1]     [1,-1]
     *
     *   - Triangle
     *        [0,1]
     *         (2)
     *          |   .
     *          |     .
     *          |        .
     *         (0) ------ (1)
     *        [0,0]      [1,0]
     *
     */
    enum class ElemType {
        Line, Square, Triangle
    };

    /***************     BASIS TYPES     ********************/
    enum class BFType {
        Lagrangian, // Polynomial basis such that Ni(xj) = delta_ij
        Nedelec     // Polynomial vector basis such that integral_ej { Ni · tj } = delta_ij
    };

}

#endif //UPA_ELEMENTDEFINITIONS_H
