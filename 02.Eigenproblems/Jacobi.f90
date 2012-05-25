
!!!    This program finds eigenvalues and eigenvectors of a real symmetric matrix by using Jacobi Method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       integer, parameter :: n=4
       integer i,j
       real, parameter:: eps=1.0e-09
       real egv
       real :: a(n,n),x(n,n)
       a(1,:)=(/1.0, 2.0, 1.0, 2.0/)
       a(2,:)=(/2.0, 2.0, -1.0, 1.0/)
       a(3,:)=(/1.0, -1.0, 1.0, 1.0/)
       a(4,:)=(/2.0, 1.0, 1.0, 1.0/)

       write(*,*) 'Jacobi Method:'
       write(*,*) 'All the elements of matrix A:'
       write(*,'(4f9.4)') ((a(i,j),j=1,n),i=1,n)
       call Jacobi(a,x,n,eps)
       do i=1,n
              write(*,*) '*************************************'
              write(*,*)
              write (*,*) 'Eigenvalues:'
              write (*,'(4f9.4)') a(i,i)
              write(*,*) 'Eigenvector:'
              write(*,15) (x(j,i),j=1,n)
15            format(1x,'[',f9.6,',',f9.6,',',f9.6,',',f9.6,']')
       enddo

       stop
       end program main



       subroutine Jacobi(a,x,n,eps)
       implicit none
       integer i, j, k, n
       real :: a(n,n),x(n,n)
       real :: eps, b2, bar
       real :: beta, coeff, c, s, cs, sc

!!!    initialize x(i,j)=0, x(i,i)=1
       x = 0.0
       do i=1,n
         x(i,i) = 1.0
       end do
!!!    find the sum of all off-diagonal elements squared
       b2 = 0.0
       do i=1,n
              do j=1,n
                     if (i.NE.j) b2 = b2 + a(i,j)**2
              end do
       end do

       if (b2 <= eps) return

!!!    average for off-diagonal elements /2
       bar = 0.5*b2/float(n*n)
       do while (b2.GT.eps)
              do i=1,n-1
                     do j=i+1,n
                            if (a(j,i)**2 .LT. bar) cycle        !!! do not touch small elements
                            b2 = b2 - 2.0*a(j,i)**2
                            bar = 0.5*b2/float(n*n)
!!!                         calculate coefficient c and s for Givens matrix
                            beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
                            coeff = 0.5*beta/sqrt(1.0+beta**2)
                            s = SQRT(MAX(0.5+coeff,0.0))
                            c = SQRT(MAX(0.5-coeff,0.0))
!!!                          recalculate rows i and j
                            do k=1,n
                                   cs =  c*a(i,k)+s*a(j,k)
                                   sc = -s*a(i,k)+c*a(j,k)
                                   a(i,k) = cs
                                   a(j,k) = sc
                            end do
!!!                         new matrix a_{k+1} from a_{k}, and eigenvectors
                            do k=1,n
                                   cs =  c*a(k,i)+s*a(k,j)
                                   sc = -s*a(k,i)+c*a(k,j)
                                   a(k,i) = cs
                                   a(k,j) = sc
                                   cs =  c*x(k,i)+s*x(k,j)
                                   sc = -s*x(k,i)+c*x(k,j)
                                   x(k,i) = cs
                                   x(k,j) = sc
                            enddo
                     enddo
              enddo
       enddo

       return
       end subroutine Jacobi
