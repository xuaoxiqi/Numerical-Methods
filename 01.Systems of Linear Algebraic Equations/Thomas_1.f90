
!!!    This program sloves tridiagonal linear equation A*x=b
!!!    Sparse matrix A is Tridiagonal Matrix
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

! Example refer to :
!
!       | 4   1   0   0   0 | |x1|   |  1.0 |
!       | 1   4   1   0   0 | |x2|   |  0.5 |
!       | 0   1   4   1   0 | |x3| = | -1.0 |
!       | 0   0   1   4   1 | |x4|   |  3.0 |
!       | 0   0   0   1   4 | |x5|   |  2.0 |
!
! Solution is :
!   x1=0.2     x2=0.2     x3=-0.5    x4=0.8     x5=0.3

    program main
    integer, parameter :: n=5
    real(8) :: A(n), B(n), C(n), F(n)
    real(8) :: X(n)

    data (A(i), i=1,n) /0.0d0, 1.0d0,  1.0d0, 1.0d0, 1.0d0/
    data (B(i), i=1,n) /4.0d0, 4.0d0,  4.0d0, 4.0d0, 4.0d0/
    data (C(i), i=1,n) /1.0d0, 1.0d0,  1.0d0, 1.0d0, 0.0d0/
    data (F(i), i=1,n) /1.0d0, 0.5d0, -1.0d0, 3.0d0, 2.0d0/
    X = 0.0d0

    write(*,*) 'Thomas Method:'
    write(*,*) "ax=b"
    write(*,*) "Coefficient for a(n*n) and b(n*1)"
    write(*,*)
    do i=1,n
        write(*,100) A(i),B(i),C(i),F(i)
    enddo
    call Thomas(A,B,C,F,X,n)

    write(*,*)
    write(*,*) 'Solution x:'
    write(*,'(f11.5)') (X(i),i=1,n)

100 format(2x,f7.3,f7.3,f7.3,f18.3)
    stop
    end program main



    subroutine Thomas(coeffA,coeffB,coeffC,coeffF,X,n)
    implicit none
    integer :: n, k
    real(8) :: coeffA(n), coeffB(n), coeffC(n), coeffF(n), X(n)
    real(8) :: A(n), B(n), C(n), F(n)
    real(8) :: t

    A = coeffA
    B = coeffB
    C = coeffC
    F = coeffF

    C(1) = C(1)/B(1)
    F(1) = F(1)/B(1)
    do k=2,n
        t = B(k)-C(k-1)*A(k)
        C(k)=C(k)/t
        F(k)=( F(k)-F(k-1)*A(k) )/t
    enddo
    do k=n-1,1,-1
        F(k)=F(k)-C(k)*F(k+1)
    enddo

    X = F

    return
    end subroutine Thomas
