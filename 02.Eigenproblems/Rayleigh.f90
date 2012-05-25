
!!!    This program finds eigenvalues and eigenvectors of a real symmetric matrix by using Rayleigh Method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       integer, parameter :: n=3
       integer i, j
       real egv
       real :: a(n,n),b(n)
       a(1,:)=(/2, 3, 2/)
       a(2,:)=(/10, 3, 4/)
       a(3,:)=(/3, 6, 1/)
       b =(/0.0, 0.0, 1.0/)

       write(*,*) 'Rayleigh Method:'
       write(*,*) 'All the elements of matrix A:'
       write(*,'(3f9.4)') ((a(i,j),j=1,n),i=1,n)
       write(*,*) 'Initial eigenvector:'
       write(*,15) b
       call Rayleigh(a,b,egv,n)
       write(*,*) '***********************************'
       write(*,*) 'Maximum eigenvalue:'
       write(*,*) egv
       write(*,*) 'Eigenvector:'
       write(*,15) b
15     format(1x,'[',f9.6,',',f9.6,',',f9.6,']')

       stop
       end program main



       subroutine Rayleigh(a,b,egv,n)
       implicit none
       integer i, j, n, itc
       real, parameter:: eps=1.0e-06
       real :: a(n,n), b(n), c(n), egv
       real :: e, e1, r0, r, fm, fz

       r0=0
       do itc=1,50
              e=0.0
              do i=1,n
                     c(i)=0.0
                     do j=1,n
                            c(i)=c(i)+a(i,j)*b(j)
                     enddo
              enddo            !!!        c=ax
              do i=1,n
                     if(ABS(c(i)).LT.e) exit
                     e=ABS(c(i))
                     e1=c(i)         !find the maximum value of c(i)
              enddo
              fm=0.0
              fz=0.0
              do i=1,n
                     fz=fz+b(i)*c(i)      !!!fz = x'ax
                     fm=fm+b(i)*b(i)      !!!fm = x'x
              enddo
              do i=1,n
                     b(i)=c(i)/e1         !!!normalization
              enddo
              r=fz/fm     !!!      Rayleigh = x'ax / x'x
              if(ABS(r-r0)/ABS(r).LT.eps) exit
              r0=r
       enddo
       egv = r

       return
       end
