

#ifndef UPA_ITERATOR_H
#define UPA_ITERATOR_H

namespace upa {

    /** General Iterator class
     *
     *  Each iterator has two main functions:
     *      1) A constructor, which sets up the internal variables.
     *      2) An iteration function, which should have no further setup required.
     */
    class Iterator {

    public:
        virtual ~Iterator() = default;
        virtual void iterate(double* x_in, double* b, double* x_out) = 0;

    };

    /** Enum Class for different iterators **/
    enum class Iterator_Type {
        Jacobi
    };

}


#endif //UPA_ITERATOR_H
