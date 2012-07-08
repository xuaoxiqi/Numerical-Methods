
!!!    This program solves Burgers Equation using 3rd-order finite volume method
!!!    Using Godunov numerical flux, with TVB Limiter
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       integer, parameter :: N=640
       real(8), parameter :: Pi=3.1415926535d0
       integer :: i, nt
       real(8) :: alpha, dx, dt, t, u_error
       real(8) :: u_avg(0:N-1), X(0:N), u_exact(0:N)


       alpha = 0.5d0
       dx = 2.0d0*Pi/N
       dt = alpha*dx
       t = 0.1d0*Pi
       nt = NINT(t/dt)

       do i=0,N
              X(i) = i*dx
       enddo


       print*, ''
       print*,'3rd FVM SCHEME(with TVB Limiter):'
       print*,'****************************************'


       call exact(N,dx,X,t,u_exact)

       call ENO(N,dx,X,dt,nt,u_avg)

       !!! calculate error
       u_error = 0.0d0
       do i=0,N-1
              u_error = ABS(u_exact(i)-u_avg(i))+u_error
       enddo
       u_error = u_error/N

       write(*,*) 'Number of points is:',N
       write(*,*) 'L1 Normal is:', u_error

       open(unit=01,file='./results.dat',status='unknown')
       write(01,101)
       write(01,102)
       write(01,103) N+1

       do i = 0,N-1
                     write(01,100) X(i), u_avg(i),u_exact(i)
       enddo
       write(01,100) X(N), u_avg(0),u_exact(N)


       close(01)
       print*,'****************************************'


100    format(2x,10(e12.6,'      '))
101    format('Title="Burgers Equation"')
102    format('Variables=x,u_num,u_exact')
103    format('zone',1x'i=',1x,i5,2x,'f=point')
       stop
       end program main



       subroutine exact(N,dx,X,t,u_exact)
       implicit none
       integer :: i, N, j
       real(8) :: dx, x0, temp, t
       real(8) :: X(0:N)
       real(8) :: u_exact(0:N), u_t(3), u_x(3), weight(3)

       u_t(1) = -0.774596669241483d0
       u_t(2) = 0.0d0
       u_t(3) = 0.774596669241483d0

       weight(1) = 0.5555555555556d0
       weight(2) = 0.8888888888889d0
       weight(3) = 0.5555555555556d0

       do i=0,N
              do j=1,3
              u_x(j) = X(i)-dx*0.5d0+dx*(1.0d0+u_t(j))/2.d0
              enddo
              u_exact(i)=0.d0
              do j=1,3
                     x0 = u_x(j)
                     do
                            temp = x0
                            x0 = x0-(x0+t*(1.0d0/3.0d0+2.0d0/3.0d0*SIN(x0))-u_x(j))/(1.0d0+t*2.0d0/3.0d0*COS(x0))
                            if(ABS(temp-x0).GT.1e-13) then
                                   cycle
                            else
                                   exit
                            endif
                     enddo
                     u_exact(i) = u_exact(i)+weight(j)*(1.0d0/3.0d0+2.0d0/3.0d0*SIN(x0))*0.5d0
              enddo

              !print*,'x0=',X(i),'u_exact=',u_exact(i)
       enddo



       return
       end subroutine exact




       subroutine ENO(N,dx,X,dt,nt,u_avg)
       implicit none
       integer :: N, i, j, nt
       real(8) :: dx, dt, ul, ur
       real(8) :: X(0:N), u_avg(0:N-1), du_avg(0:N), f(N), u1(0:N-1), u2(0:N-1)

       do i=1,N-1
              u_avg(i) = 1.0d0/3.0d0+4.0d0/3.0d0/dx*SIN(X(i))*SIN(0.5d0*dx)
       enddo

       u_avg(0) = 1.0/3.0

       do j=1,nt

              call Godunov(u_avg,du_avg,dx,N)
              do i=0,N-1
                     u1(i) = u_avg(i)+dt*du_avg(i)
              enddo

              call Godunov(u1,du_avg,dx,N)
              do i=0,N-1
                     u2(i) = 3.0d0/4.0d0*u_avg(i)+1.0d0/4.0d0*(u1(i)+dt*du_avg(i))
              enddo

              call Godunov(u2,du_avg,dx,N)
              do i=0,N-1
                     u_avg(i) = 1.0d0/3.0d0*u_avg(i)+2.0d0/3.0d0*(u2(i)+dt*du_avg(i))
              enddo
       enddo


!       do i=0,N-1
!              print*,'x0=',X(i),'u(i)=',u_avg(i)
!       enddo
!       print*,'x0=',X(N),'u(i)=',u_avg(0)

       return
       end subroutine ENO


       subroutine Godunov(u_avg, du_avg, dx, N)
       implicit none
       integer :: i, N
       real(8) :: u_avg(0:N-1), du_avg(0:N-1), ul, ur, f(N), dx

       do i=0,N-1
              ul = -1.0d0/6.0d0*u_avg(MOD(i-1+N,N))+5.0/6.0*u_avg(i)+1.0d0/3.0d0*u_avg(MOD(i+1,N))
              ur = 1.0d0/3.0d0*u_avg(i)+5.0d0/6.0d0*u_avg(MOD(i+1,N))-1.0d0/6.0d0*u_avg(MOD(i+2,N))

              call TVB_Limiter(ul,ur,u_avg(MOD(i-1+N,N)),u_avg(i),u_avg(MOD(i+1,N)),u_avg(MOD(i+2,N)),dx)

              if(ul.LE.ur) then
                     f(i+1) = MIN(0.5d0*ul*ul,0.5d0*ur*ur)
                     if(ul*ur.LT.0.0d0) f(i+1) = 0.0d0
              else
                     f(i+1) = MAX(0.5d0*ul*ul,0.5d0*ur*ur)
              endif
       enddo

       do i=1,N-1
              du_avg(i) = -(f(i+1)-f(i))/dx
       enddo
       du_avg(0) = -(f(1)-f(N))/dx
       end subroutine Godunov



       subroutine TVB_Limiter(ul,ur,u_avg_j0,u_avg_j1,u_avg_j2,u_avg_j3,dx)
       implicit none
       real(8) :: ul, ur, u_avg_j0, u_avg_j1, u_avg_j2, u_avg_j3, u_tilde, u_dtilde, dx, M

       M = 5.0d0

       u_tilde = ul-u_avg_j1
       u_dtilde = u_avg_j2-ur

       if (ABS(u_tilde).GT.M*dx*dx) then
              if ((SIGN(1.0d0,u_tilde).EQ.SIGN(1.0d0,u_avg_j2-u_avg_j1)).AND. &
                 &(SIGN(1.0d0,u_tilde).EQ.SIGN(1.0d0,u_avg_j1-u_avg_j0))) then
                     u_tilde = MIN(ABS(u_tilde),ABS(u_avg_j2-u_avg_j1),ABS(u_avg_j1-u_avg_j0))*SIGN(1.0d0,u_tilde)
              else
                     u_tilde = 0.0d0
              endif
       endif

       if (ABS(u_dtilde).GT.M*dx*dx) then
              if ((SIGN(1.0d0,u_dtilde).EQ.SIGN(1.0d0,u_avg_j3-u_avg_j2)).AND. &
                 &(SIGN(1.0d0,u_tilde).EQ.SIGN(1.0d0,u_avg_j2-u_avg_j1))) then
                     u_dtilde = MIN(ABS(u_dtilde),ABS(u_avg_j3-u_avg_j2),ABS(u_avg_j2-u_avg_j1))*SIGN(1.0d0,u_dtilde)
              else
                     u_dtilde = 0.0d0
              endif
       endif

       ul = u_avg_j1+u_tilde
       ur = u_avg_j2+u_dtilde

       return
       end subroutine TVB_Limiter



