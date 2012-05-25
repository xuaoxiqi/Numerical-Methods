
!!!    This program finds the maximum eigenvalue of a matrix and its eigenvector for a given initial eigenvector by using power method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       integer, parameter :: n=3
       integer i,j
       real egv
       real :: a(n,n), b(n), c(n)
       a(1,:)=(/2, 3, 2/)
       a(2,:)=(/10, 3, 4/)
       a(3,:)=(/3, 6, 1/)
       b = (/0.0, 0.0, 1.0/)

       write(*,*) 'Power Method:'
       write(*,*) 'All the elements of matrix A:'
       write(*,'(3f9.4)') ((a(i,j),j=1,n),i=1,n)
       write(*,*) 'Initial eigenvector:'
       write(*,15) b
       call power(a,b,c,egv,n)
       write(*,*) '***********************************'
       write(*,*) 'Maximum eigenvalue:'
       write(*,*) egv
       write(*,*) 'Eigenvector:'
       write(*,15) b
15     format(1x,'[',f9.6,',',f9.6,',',f9.6,']')


       stop
       end program main



       subroutine power(a,b,c,egv,n)
       implicit none
       integer n,itc,i,j
       real :: a(n,n),b(n),c(n)
       real :: e,e1,egv

       do itc=1,50
              e=0.0
              do i=1,n
                     c(i)=0.0
                     do j=1,n
                            c(i)=c(i)+a(i,j)*b(j)
                     end do
              end do
              do i=1,n
                     if(ABS(c(i)).LT.e) cycle
                     e=ABS(c(i))
                     e1=c(i)
              end do
              do i=1,n
                     b(i)=c(i)/e1
              end do
       end do
       egv=e1

       return
       end subroutine power
