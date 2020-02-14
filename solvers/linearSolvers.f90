
MODULE linearSolvers
implicit NONE
contains

! #####################################################
! ##################      CHOLESKY     ################
! #####################################################

! Returns the Cholesky decomposition of S in D
! To have a Cholesky decomposition, S must be a 
! square symmetric and positive definite matrix
SUBROUTINE Chol(n , S , D)
    implicit NONE
    INTEGER :: n ! Size of the matrixes
    INTEGER :: i , j , k ! iterators
    REAL :: S(n , n) , D(n , n)
    REAL :: sum

    ! We fill the lower diagonal matrix
    DO i = 1 , n
        DO j = 1 , i
            sum = 0.0
            IF (j == i) THEN ! Diagonal elements
                DO k = 1, j-1
                    sum = sum + D(j,k)*D(j,k)
                END DO
                D(j,j) = SQRT(S(j,j)-sum)
            ELSE ! Non-Diagonal elements
                DO k = 1, j-1
                    sum = sum + D(i,k)*D(j,k)
                END DO
                D(i,j) = (S(i,j)-sum)/D(j,j)
            END IF 
        END DO 
    END DO
END SUBROUTINE Chol


SUBROUTINE solveChol(n , L , b , x)
    implicit NONE
    INTEGER n
    REAL L(:,:), b(:), x(:)
    REAL sum, y(n)

    ! First, we solve L·y = b (Lower triangular system)
    DO i = 1 , n
        sum = 0.0
        DO j = 1 , i-1
            sum = sum + L(i,j)*y(j)
        END DO
        y(i) = (b(i) - sum) / L(i,i)
    END DO

    ! Finally, we solve L'·x = y (Upper triangular system)
    DO i = n , -1 , 1
        sum = 0.0
        DO j = i+1 , n
            sum = sum + L(j,i)*x(j)
        END DO
        x(i) = (y(i) - sum) / L(i,i)
    END DO

END SUBROUTINE solveChol



! #####################################################
! ###########      LU WITH PARTIAL PIVOTING     #######
! #####################################################

! PA = LU with row pivoting gauss decomposition of A (square nxn matrix)
! L is a unity lower triangular matrix, Lii = 1.0 implicitly
! U is an upper triangular matrix
SUBROUTINE LU( n , A , P)
    implicit NONE
    INTEGER :: n
    REAL :: A(n , n) , L(n , n) , U(n , n) 
    INTEGER :: P(n)

    INTEGER :: i , j , k, pivot
    REAL :: auxElem , auxIndex

    ! Initialise matrixes
    DO i = 1, n
        P(i) = i
    END DO

    ! Main Loop
    DO k = 1 , n ! Loop over columns

        ! Find maximum of row
        max = abs(A(k , k))
        pivot = k
        DO j = k , n
            IF (max < abs(A(k , j))) THEN
                pivot = j
                max = abs(A(k , j))
            END IF
        END DO

        IF (max < 1.e-3) THEN ! the hole row is null
            stop ! Matrix is singular
        END IF

        IF (k .NE. pivot) THEN ! Do we need to swap rows?
            ! swap rows
            DO j = 1 , n
                auxElem = A(k , j)
                A(k , j) = A(pivot , j)
                A(pivot , j) = auxElem
            END DO
            ! write change in permutation vector
            auxIndex = P(i)
            P(i) = P(maxj)
            P(maxj) = auxIndex
        END IF

        ! We calculate the kth row of U
        ! Ukj = Akj - sum_{i=1,...k-1} Uij Lki  , j >= k
        DO j = k , n
            sum = 0.0
            DO i = 1 , k-1
                sum = sum + A(i,j)*A(k,i)
            END DO
            A(k,j) = A(k,j) - sum
        END DO

        ! We calculate the kth column of L
        ! Lik = ( Aik - sum_{j=1,...k-1} Ujk Lij ) / Ujj , k < i
        ! Lkk = 1.0 implicitly
        DO i = k+1 , n
            sum = 0.0
            DO j = 1, k-1
                sum = sum + A(j,k)*A(i,j)
            END DO
            A(i,k) = (A(i,k) - sum) / A(k,k)
        END DO

    END DO ! End Main Loop
END SUBROUTINE LU


SUBROUTINE solveLU()
    implicit NONE

END SUBROUTINE solveLU



! #####################################################
! ###########      QR WITH PARTIAL PIVOTING     #######
! #####################################################

! Decomposes AP = QR using the Householder method with column pivoting
! A is of size m*n , with m >= n
! A gets overwritten during the procedure
SUBROUTINE QR( m , n , A , betaVector , P)
    implicit NONE
    INTEGER :: m , n
    REAL :: A(: , :) , betaVector(:)
    INTEGER :: P(:)

    REAL , allocatable :: v(:) , x(:)
    REAL :: beta , colnorms(n) , swapReal
    INTEGER :: i , j , k , pivot , swapInt

    ! Setup for pivoting mechanism 
    DO j = 1 , n
        P(j) = j
        colnorms(j) = 0.0
        DO i = 1 , m
            colnorms(j) = colnorms(j) + A(i , j)*A(i , j)
        END DO
    END DO

    ! Main loop, iterating through columns
    DO j = 1 , n

        ! Calculate pivot
        pivot = j
        DO k = j+1 , n
            IF (colnorms(k) > colnorms(pivot)) pivot = k
        END DO

        IF (colnorms(pivot) < 1.e-14) stop
        IF (pivot .NE. j) THEN ! swap columns j and pivot
            swapInt = P(pivot)
            P(pivot) = P(j)
            P(j) = swapInt

            swapReal = colnorms(pivot)
            colnorms(pivot) = colnorms(j)
            colnorms(j) = swapReal

            DO k = 1 , m
                swapReal = A(k , pivot)
                A(k , pivot) = A(k , j)
                A(k , j) = swapReal
            END DO
        END IF
        ! End of pivoting: column j is the one to modify

        ! Calculate Householder matrix of column A(j:m , j)
        IF (allocated(x)) deallocate(x)
        IF (allocated(v)) deallocate(v)
        allocate(x(m - j + 1))
        allocate(v(m - j + 1))
        DO i = 1 , m-j+1
            x(i) = A(j+i-1 , j)
        END DO
        CALL Householder(m-j+1 , x , v , beta) ! H = I - beta * v·v'

        ! A(j:m , j:n) = A(j:m , j:n) - beta*v·v'·A(j:m , j:n)
        ! We will use x as buffer for v'·A(j:m , j:n)
        deallocate(x)
        allocate(x(n - j + 1))
        ! Calculate x = v'·A(j:m , j:n)
        DO k = 1 , n-j+1 
            x(k) = 0.0
            DO i = 1 , m-j+1
                x(k) = x(k) + v(i)*A(j+i-1 , j+k-1)
            END DO
        END DO
        ! Modify A(j:m , j:n)
        DO k = 1 , n-j+1
            DO i = 1 , m-j+1
                A(j+i-1 , j+k-1) = A(j+i-1 , j+k-1) - beta*v(i)*x(k)
            END DO
        END DO

        ! Copy v into A to save it
        DO i = 2 , m-j+1 ! We do not save v(1) since it is always 1
            A(j+i-1 , j) = v(i)
        END DO
        ! Save beta
        betaVector(j) = beta

        ! Update colnorm
        DO k = j+1 , n
            colnorms(k) = colnorms(k)-A(j , k)*A(j , k)
        END DO

    END DO
END SUBROUTINE QR


SUBROUTINE solveQR()
    implicit NONE

END SUBROUTINE solveQR

! Computes the householder matrix of x
! Input: n , x (x has size n)
! Output: v , beta such that 
!          - v(1) = 1 , v has size n
!          - beta is a scalar
!          - P = I - beta*v·v' is orthogonal
!          - P·x = norm(x)*e_1
SUBROUTINE Householder(n , x , v , beta)
    implicit NONE
    INTEGER :: n
    REAL, allocatable :: x(:) , v(:) 
    REAL :: beta

    REAL :: s , mu
    INTEGER :: i , j

    v(1) = 1.0
    s = 0.0
    DO i = 2 , n
        s = s + x(i)*x(i)
        v(i) = x(i)
    END DO

    IF (s < 1.e-14) THEN
        beta = 0.0
    ELSE 
        mu = sqrt(x(1)*x(1) + s)
        IF (x(1) .LE. 0) THEN
            v(1) = x(1) - mu
        ELSE 
            v(1) = -s/(x(1) + mu)
        END IF
        beta = 2.0*v(1)*v(1)/(s + v(1)*v(1))
        DO i = 2 , n
            v(i) = v(i)/v(1)
        END DO
        v(1) = 1.0
    END IF
END SUBROUTINE Householder

END MODULE
