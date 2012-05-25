
!!!    This program sloves tridiagonal linear equation A*x=b
!!!    Diagonally dominant matrix(sparse matrix) coefficient is stored in n*3 matrix
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none

       integer, parameter :: n=7  !!! n is number of equations
       integer :: i,j
       real :: a(n,3), b(n), x(n)  !!! x(i) is solution vector
       data (a(1,j), j=1,3) / 0.0, -2.25, 1.0/  !!! a is coefficient matrix A(i,3)
       data (a(2,j), j=1,3) / 1.0, -2.25, 1.0/
       data (a(3,j), j=1,3) / 1.0, -2.25, 1.0/
       data (a(4,j), j=1,3) / 1.0, -2.25, 1.0/
       data (a(5,j), j=1,3) / 1.0, -2.25, 1.0/
       data (a(6,j), j=1,3) / 1.0, -2.25, 1.0/
       data (a(7,j), j=1,3) / 1.0, -2.25, 0.0/
       data (b(i),i=1,n) / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -100.0 /  !!! b(i) is right hand side vector

       write(*,*) 'Thomas Method:'
       write(*,*) 'All the elements of matrix A:'
       write(*,'(3f11.5)') ((a(i,j),j=1,3),i=1,n)
       write(*,*)
       write(*,*) 'Vector b:'
       write(*,'(f11.5)') (b(i),i=1,n)
       call Thomas(a,b,x,n)
       write(*,*)
       write(*,*) 'Solution x:'
       write(*,'(f11.5)') (x(i),i=1,n)

       stop
       end program main



       subroutine Thomas(a,b,x,n)
       implicit none

       integer :: n, i
       real :: em
       real :: a(n,3), b(n), x(n)
!!!    forward elimination
       do i=2,n
              em=a(i,1)/a(i-1,2)
              a(i,1)=em
              a(i,2)=a(i,2)-em*a(i-1,3)
              b(i)=b(i)-a(i,1)*b(i-1)
       enddo
!!!    back substitution
       x(n)=b(n)/a(n,2)
       do i=n-1,1,-1
              x(i)=(b(i)-a(i,3)*x(i+1))/a(i,2)
       enddo

       return
       end subroutine Thomas
