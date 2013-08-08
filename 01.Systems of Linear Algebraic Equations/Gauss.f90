
!!!    This program solves a system of linear equation A*x=b by using Gauss elimination(with scaled pivoting)
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

    program main
    implicit none
    integer, parameter :: N=4
    integer, parameter :: dp=kind(0.d0)
    real(dp) :: a(N,N), b(N), x(N)
    integer i,j
    ! matrix A
    data (a(1,i), i=1,N) /  2.0_dp, 0.0_dp, 0.0_dp, 6.0_dp /
    data (a(2,i), i=1,N) /  3.0_dp, 7.0_dp, 0.0_dp, 0.0_dp /
    data (a(3,i), i=1,N) /  0.0_dp, 0.0_dp, 1.5_dp, 0.0_dp /
    data (a(4,i), i=1,N) /  0.0_dp, 10.0_dp, 0.0_dp, 8.0_dp /
    !matrix b
    data (b(i),   i=1,N) / 13.0_dp, 8.5_dp, 2.25_dp, 26.0_dp /

    write(*,*) 'Gauss elimination with scaling and pivoting method:'
    write(*,*) 'All the elements of matrix A:'
    write(*,'(4f9.4)') ((a(i,j),j=1,N),i=1,N)
    write(*,*)
    write(*,*) 'Vector b:'
    write(*,'(f9.4)') (b(i),i=1,N)

    call Gauss(a,b,x,N)

    write(*,*)
    write(*,*) 'Solution x:'
    write(*,'(f9.4)') (x(i),i=1,N)

    stop
    end program main



    subroutine Gauss(a,b,x,N)
    implicit none
    integer N
    integer, parameter :: dp=kind(0.d0)
    real(dp) :: a(N,N), b(N), x(N)
    real(dp) :: s(N)
    real(dp) :: c, pivot, store
    integer i, j, k, l
!!!    begin forward elimination
    do k=1, N-1
!!!     scaling
!!!     s(i) will have the largest element from row i
        do i=k,N                       !!! loop over rows
            s(i) = 0.0_dp
            do j=k,N                    !!! loop over elements of row i
                s(i) = MAX(s(i),ABS(a(i,j)))
            end do
        end do
!!!           pivoting 1
!!!           find a row with the largest pivoting element
        pivot = ABS(a(k,k)/s(k))
        l = k
        do j=k+1,N
            if(ABS(a(j,k)/s(j)).GT.pivot) then
                pivot = abs(a(j,k)/s(j))
                l = j
            end if
        end do

!!!           Check if the system has a sigular matrix
        if(pivot.EQ.0.0_dp) then
            write(*,*) ' The matrix is sigular '
            return
        end if

!!!    pivoting 2 interchange rows k and l (if needed)
        if (l /= k) then
            do j=k,n
                store = a(k,j)
                a(k,j) = a(l,j)
                a(l,j) = store
            end do
            store = b(k)
            b(k) = b(l)
            b(l) = store
        end if

!!!    the elimination (after scaling and pivoting)
        do i=k+1,n
            c = a(i,k)/a(k,k)
            a(i,k) = 0.0_dp
            b(i) = b(i)-c*b(k)
            do j=k+1,n
                a(i,j) = a(i,j)-c*a(k,j)
            enddo
        enddo
    enddo
!!!    back substiturion
    x(n) = b(n)/a(n,n)
    do i=n-1,1,-1
        c = 0.0_dp
        do j=i+1,n
            c = c+a(i,j)*x(j)
        enddo
        x(i) = (b(i)- c)/a(i,i)
    enddo

    return
    end subroutine Gauss
