

! Test for the linear solvers
PROGRAM testSolvers

    write(*,*) 'Testing Cholesky : '
    call testChol

    write(*,*) 'Testing LU : '
    call testLU

END PROGRAM testSolvers


SUBROUTINE testChol()
    USE linearSolvers
    implicit NONE

    INTEGER n, i
    REAL , allocatable :: A(:,:) , b(:) , x(:) , L(:,:)

    ! Inputs
    n = 3
    allocate(A(n,n) , b(n) , L(n,n) , x(n))

    A(1,:) = (/  2.0 , -1.0 ,  0.0 /)
    A(2,:) = (/ -1.0 ,  2.0 , -1.0 /)
    A(3,:) = (/  0.0 , -1.0 ,  2.0 /)

    b = (/ 1 , 2 , 3/)

    ! Solve using LU
    call Chol(n, A, L)
    call solveChol(n , L , b , x)

    ! Output
    write(*,*) 'n = ' , n
    write(*,*) 'x = ' , x
    write(*,*) 'L = '
    DO i = 1 , n
        write(*,*) L(i,:)
    END DO

    deallocate(A, b, L, x)

END SUBROUTINE testChol


SUBROUTINE testLU()
    USE linearSolvers
    implicit NONE

    INTEGER n, i
    REAL , allocatable :: A(:,:) , b(:) , x(:)
    INTEGER , allocatable :: P(:)

    ! Inputs
    n = 3
    allocate(A(n,n) , b(n) , P(n) , x(n))

    A(1,:) = (/  2.0 , -1.0 ,  0.0 /)
    A(2,:) = (/ -1.0 ,  2.0 , -1.0 /)
    A(3,:) = (/  0.0 , -1.0 ,  2.0 /)

    b = (/ 1 , 2 , 3/)

    ! Solve using LU
    call LU(n, A, P)
    call solveLU(n , A , P , b , x)

    ! Output
    write(*,*) 'n = ' , n
    write(*,*) 'x = ' , x
    write(*,*) 'LU = '
    DO i = 1 , n
        write(*,*) A(i,:)
    END DO
    write(*,*) 'P = ', P

    deallocate(A, b, P, x)

END SUBROUTINE testLU

