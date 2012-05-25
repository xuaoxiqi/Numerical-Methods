
!!!    This program finds eigenvalues of a matrix and its eigenvector for a given initial eigenvalue by using inverse power method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       integer, parameter :: n=4
       integer i,j
       real egv
       real :: a(n,n),x(n)
       real :: y(n)=(/9.7, 0.6, -0.6, -1.7/)
       a(1,:)=(/3.0, 2.0, 3.0, 4.0/)
       a(2,:)=(/2.0, 1.0, 1.0, 1.0/)
       a(3,:)=(/3.0, 1.0, 2.0, 3.0/)
       a(4,:)=(/4.0, 1.0, 3.0, 2.0/)

       write(*,*) 'Inverse power Method:'
       write(*,*) 'All the elements of matrix A:'
       write(*,'(3f9.4)') ((a(i,j),j=1,n),i=1,n)
       write(*,*) 'Initial eigenvalue:'
       write(*,'(f9.4)') y

       do i=1,n
              write(*,*) '*************************************'
              write(*,*)
              egv=y(i)
              call ipmcv(a,x,egv,n)
              write(*,*) 'Eigenvalue:',egv
              write(*,*) 'Eigenvector:'
              write(*,15) x
15            format(1x,'[',f9.6,',',f9.6,',',f9.6,',',f9.6,']')
              write(*,*)
       enddo

       stop
       end program main

       subroutine ipmcv(a,x,egv,n)
       implicit none
       integer n, i, j, k, itc
       real :: e, e0, e1, egv
       real :: a(n,n), b(n,n), x(n)

       b=a
       do i=1,n
              b(i,i)=b(i,i)-egv
       end do
       do j=1,n-1
              do k=j+1,n
                     b(j,k)=b(j,k)/b(j,j)
                     do i=j+1,n
                            b(i,k)=b(i,k)-b(i,j)*b(j,k)
                     end do
              end do
       end do
       do i=1,n
              x(i)=1.0
       end do
       do itc=1,50
              do i=n-1,1,-1
                     do j=n,i+1,-1
                            x(i)=x(i)-b(i,j)*x(j)
                     enddo
              enddo
              e=0.0
              do i=1,n
                     if(ABS(x(i)).LT.e) cycle
                     e=ABS(x(i))
                     e1=x(i)
              enddo
              do i=1,n
                     x(i)=x(i)/e1
              enddo
              if(ABS(e-e0).LT.5.0E-7*ABS(e)) exit
              e0=e
              x(1)=x(1)/b(1,1)
              do i=2,n
                     do j=1,i-1
                            x(i)=x(i)-b(i,j)*x(j)
                     end do
                     x(i)=x(i)/b(i,i)
              enddo
       enddo
       egv=egv+1.0/e1

       return
       end subroutine ipmcv
