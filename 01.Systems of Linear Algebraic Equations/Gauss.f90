
!!!    This program solves a system of linear equation A*x=b by using Gauss elimination(with scaled pivoting)
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       integer, parameter :: n=4
       double precision a(n,n), b(n), x(n)
       integer i,j
       ! matrix A
       data (a(1,i), i=1,4) /  2.0, 0.0, 0.0, 6.0 /
       data (a(2,i), i=1,4) /  3.0, 7.0, 0.0, 0.0 /
       data (a(3,i), i=1,4) /  0.0, 0.0, 1.5, 0.0 /
       data (a(4,i), i=1,4) /  0.0, 10.0, 0.0, 8.0 /
       ! matrix b
       data (b(i),   i=1,4) / 13.0, 8.5, 2.25, 26.0 /

       write(*,*) 'Gauss elimination with scaling and pivoting method:'
       write(*,*) 'All the elements of matrix A:'
       write(*,'(4f9.4)') ((a(i,j),j=1,n),i=1,n)
       write(*,*)
       write(*,*) 'Vector b:'
       write(*,'(f9.4)') (b(i),i=1,n)
       call Gauss(a,b,x,n)
       write(*,*)
       write(*,*) 'Solution x:'
       write(*,'(f9.4)') (x(i),i=1,n)

       stop
       end program main



       subroutine Gauss(a,b,x,n)
       implicit none
       integer n
       double precision a(n,n), b(n), x(n)
       double precision s(n)
       double precision c, pivot, store
       integer i, j, k, l

!!!    begin forward elimination
       do k=1, n-1
!!!           scaling
       !!!    s(i) will have the largest element from row i
              do i=k,n                       !!! loop over rows
                     s(i) = 0.0
                     do j=k,n                    !!! loop over elements of row i
                            s(i) = MAX(s(i),ABS(a(i,j)))
                     end do
              end do
!!!           pivoting 1
!!!           find a row with the largest pivoting element
              pivot = ABS(a(k,k)/s(k))
              l = k
              do j=k+1,n
                     if(ABS(a(j,k)/s(j)) > pivot) then
                            pivot = abs(a(j,k)/s(j))
                            l = j
                     end if
              end do

!!!           Check if the system has a sigular matrix
              if(pivot.EQ.0.0) then
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
                     c=a(i,k)/a(k,k)
                     a(i,k) = 0.0
                     b(i)=b(i)- c*b(k)
                     do j=k+1,n
                            a(i,j) = a(i,j)-c*a(k,j)
                     enddo
              enddo
       enddo
!!!    back substiturion
       x(n) = b(n)/a(n,n)
       do i=n-1,1,-1
              c=0.0
              do j=i+1,n
                     c= c + a(i,j)*x(j)
              enddo
              x(i) = (b(i)- c)/a(i,i)
       enddo

       return
       end subroutine Gauss
