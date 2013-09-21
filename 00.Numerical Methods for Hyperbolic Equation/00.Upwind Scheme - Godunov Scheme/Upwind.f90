
!!!    This program solves Burgers Equation using 1st-order Upwind Scheme.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       integer, parameter :: N=640
       real(8), parameter :: Pi=3.14159265358979d0
       integer :: i, nt
       real(8) :: alpha, dx, dt, t, u_error
       real(8) :: u(0:N), X(0:N), u_exact(0:N)

       alpha = 0.5d0
       dx = Pi/DBLE(N)
       dt = alpha*dx
       t = 0.1d0*Pi
       nt = NINT(t/dt)

       do i=0,N
              X(i) = i*dx
       enddo

       write(*,*) ''
       write(*,*) 'UPWIND SCHEME:'
       write(*,*) '****************************************'

       call exact(N,X,t,u_exact)

       call FTBS(nt,N,alpha,u,X)

!!! calculate error
       u_error = 0.0d0
       do i=0,N-1
              u_error = ABS(u_exact(i)-u(i))+u_error
       enddo
       u_error = u_error/DBLE(N)

       write(*,*) 'Number of points is:',N
       write(*,*) 'L1 Normal is:', u_error

       open(unit=01,file='./result.dat',status='unknown')
       write(01,101)
       write(01,102)
       write(01,103) N+1

       do i = 0,N
                     write(01,100) X(i), u(i),u_exact(i)
       enddo

       close(01)
       print*,'****************************************'


100    format(2x,10(e12.6,'      '))
101    format('Title="Burgers Equation"')
102    format('Variables=x,u_num,u_exact')
103    format('zone',1x'i=',1x,i5,2x,'f=point')
       stop
       end program main



       subroutine exact(N,X,t,u_exact)
       implicit none
       integer :: i, N
       real(8) :: x0, temp, t
       real(8) :: X(0:N), u_exact(0:N)

       do i=0,N
              x0 = X(i)
                     do
                            temp = x0
                            x0 = x0-(x0+t*SIN(x0)-X(i))/(1.0d0+t*COS(x0))
                            if(ABS(temp-x0).GT.1e-5) then
                                   cycle
                            else
                                   exit
                            endif
                     enddo
              u_exact(i) = SIN(x0)
              !write(*,*) 'x0=',X(i),'u_exact=',u_exact(i)
       enddo

       return
       end subroutine exact



       subroutine FTBS(nt,N,alpha,u,X)
       implicit none
       integer :: i, j, nt, N
       real(8) :: alpha
       real(8) :: un(0:N), u(0:N), X(0:N)

       do i=0,N
              un(i) = SIN(X(i))
       enddo

       do j=1,nt
              u(0) = un(0)-alpha*un(0)*(un(0)-un(N-1))
              do i=1,N
                     u(i) = un(i)-alpha*un(i)*(un(i)-un(i-1))
              enddo
              un = u
       enddo

!       do i=1,N
!              write(*,*) 'x=',X(i),'u=',u(i)
!       enddo

       return
       end subroutine FTBS
