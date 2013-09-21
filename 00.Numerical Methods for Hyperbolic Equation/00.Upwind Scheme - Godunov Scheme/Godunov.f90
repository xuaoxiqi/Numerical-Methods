
!!!    This program solves Burgers Equation using Godunov Scheme.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       integer, parameter :: N=21
       real, parameter :: Pi=3.1415926535
       integer :: i, nt
       real :: alpha, dx, dt, t, u_error, error_0
       real :: u(N), X(N), u_exact(N)

       alpha = 0.5
       dx = 2*Pi/(N-1)
       dt = alpha*dx
       t = 0.1*Pi
       nt = NINT(t/dt)

       do i=1,N
              X(i) = (i-1)*dx
       enddo

       write(*,*) ''
       write(*,*) 'GODUNOV SCHEME:'
       write(*,*) '****************************************'

       call exact(N,X,t,u_exact)

       call Godunov(nt,alpha,N,X,u)

       !!!calculate error
       u_error = 0.0
       do i=1,N
              u_error = ABS(u(i)-u_exact(i))
              if(u_error.GT.error_0) error_0 = u_error
       enddo

       write(*,*) 'Number of points is:',N
       write(*,*) 'Infinity Normal is:', u_error

       open(unit=01,file='./result.dat',status='unknown')
       write(01,101)
       write(01,102)
       write(01,103) N

       do i = 1,N
                     write(01,100) X(i), u(i),u_exact(i)
       enddo

       close(01)
       print*,'****************************************'


100    format(2x,10(e12.6,'      '))
101    format('Title="Burgers Equation(Godunov Scheme)"')
102    format('Variables=x,u_num,u_exact')
103    format('zone',1x'i=',1x,i5,2x,'f=point')
       stop
       end program main



       subroutine exact(N,X,t,u_exact)
       implicit none
       integer :: i, N
       real :: x0, temp, t
       real :: X(N), u_exact(N)

       do i=1,N
              x0 = X(i)
                     do
                            temp = x0
                            x0 = x0-(x0+t*(1.0/3.0+2.0/3.0*SIN(x0))-X(i))/(1.0+t*2.0/3.0*COS(x0))
                            if(ABS(temp-x0).GT.1e-5) then
                                   cycle
                            else
                                   exit
                            endif
                     enddo
              u_exact(i) = 1.0/3.0+2.0/3.0*SIN(x0)
              !write(*,*) 'x0=',X(i),'u_exact=',u_exact(i)
       enddo

       return
       end subroutine exact


       subroutine Godunov(nt,alpha,N,X,u)
       implicit none
       integer :: N, i, j, nt
       real :: alpha
       real :: un(N), u(N), f(0:N), X(N)

       do i=1,N
              un(i) = 1.0/3.0+2.0/3.0*SIN(X(i))
       enddo

       do j=1,nt
              do i=1,N-1
                     if(un(i).LE.un(i+1)) then
                            f(i) = MIN(0.5*un(i)*un(i),0.5*un(i+1)*un(i+1))
                            if(un(i)*un(i+1).LT.0.0) f(i) = 0.0
                     else
                            f(i) = MAX(0.5*un(i)*un(i),0.5*un(i+1)*un(i+1))
                     endif
              enddo

              f(0) = f(N-1)
              f(N) = f(1)

              do i =1,N
                     u(i) = un(i)-alpha*(f(i)-f(i-1))
              enddo

              un = u
       enddo

!       do i=1,N
!              write(*,*) 'x=',X(i),'u=',u(i)
!       enddo

       return
       end subroutine Godunov
