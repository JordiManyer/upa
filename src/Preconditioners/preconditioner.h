
#ifndef UPA_PRECONDITIONER_H
#define UPA_PRECONDITIONER_H


namespace upa {

    class Preconditioner {

    public:
        virtual ~Preconditioner() = default;
        virtual void apply(double* r_in, double* r_out) = 0;

    };

    /** Enum Class for different preconditioners **/
    enum class Preconditioner_Type {
        Diagonal
    };

}

#endif //UPA_PRECONDITIONER_H
