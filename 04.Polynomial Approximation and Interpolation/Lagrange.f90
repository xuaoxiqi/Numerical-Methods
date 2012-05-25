
!!!    This program computes the interpolation by Lagrange interpolation.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       integer, parameter :: n=5
       real :: x,y
       real :: x0(n+1),y0(n+1)
       data x0/0.1, 0.2, 0.3, 0.4, 0.5, 0.6/
       data y0/0.0998, 0.1987, 0.2955, 0.3894, 0.4794, 0.5646/

       write(*,*) 'Lagrange interpolation:'
       write(*,*) 'Function f is defined by a set of n+1 points'
       write(*,*) 'x 0.1    0.2    0.3    0.4    0.5    0.6'
       write(*,*) 'y 0.0998 0.1987 0.2955 0.3894 0.4794 0.5646'
       write(*,*) 'Input x:'
       read(*,*) x
       call Lagrange(n,x0,y0,x,y)
       write(*,*) 'x=',x,' y=',y

       stop
       end program main



       subroutine Lagrange(n,x0,y0,x,y)
       implicit none
       integer :: k,j,n
       real :: p,y,x
       real :: x0(n+2),y0(n+2)

       y=0.0
       do k=0,n
              p=1.0
              do j=0,n
                     if(j.NE.k) p=p*(x-x0(j+1))/(x0(k+1)-x0(j+1))
              enddo
              y=y+p*y0(k+1)
       enddo
       end subroutine Lagrange


