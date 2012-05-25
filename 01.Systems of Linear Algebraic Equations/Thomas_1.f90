
!!!    This program sloves tridiagonal linear equation A*x=b
!!!    Sparse matrix A is Tridiagonal Matrix
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       integer, parameter :: n=5
       real :: A(n), B(n), C(n),F(n)

       data (A(i), i=1,n) /0.0, 1.0, 1.0, 1.0, 1.0/
       data (B(i), i=1,n) /4.0, 4.0, 4.0, 4.0, 4.0/
       data (C(i), i=1,n) /1.0, 1.0, 1.0, 1.0, 0.0/
       data (F(i), i=1,n) /1.0, 0.5, -1.0, 3.0, 2.0/

       write(*,*) 'Thomas Method:'
       write(*,*) 'Coefficient:'
       do i=1,n
              write(*,*) A(i),B(i),C(i),F(i)
       enddo
       call Thomas(A,B,C,F,n)
       write(*,*)
       write(*,*) 'Solution x:'
       write(*,'(f11.5)') (F(i),i=1,n)

       stop
       end program main



       subroutine Thomas(A,B,C,X,n)
       implicit none
       integer :: n, i, k
       real :: A(n), B(n), C(n), X(n)
       real :: t


       C(1) = C(1)/B(1)
       X(1) = X(1)/B(1)
       do k = 2,n
              t = B(k)-C(k-1)*A(k)
              C(k)=C(k)/t
              X(k)=( X(k)-X(k-1)*A(k) )/t
       enddo
       do k=n-1,1,-1
              X(k)=X(k)-C(k)*X(k+1)
       enddo

       return
       end subroutine Thomas
