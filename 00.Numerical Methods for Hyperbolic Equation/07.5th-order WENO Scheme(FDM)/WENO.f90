

!!!    This program solves Burgers Equation using WENO Scheme with Finite Difference Method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       integer, parameter :: N=640
       real(8), parameter :: Pi=3.1415926535897932385d0
       integer :: i, nt
       real(8) :: alpha, dx, dt, t, u_error, error_0
       real(8) :: u_avg(0:N-1)
       real(8) :: u_exact(0:N), X(0:N)

       dx = 2.0d0*Pi/DBLE(N)
       alpha = 0.01d0
       dt = alpha*dx
       t = 0.1d0*Pi
       nt = NINT(t/dt)

       do i=0,N
              X(i) = DBLE(i)*dx
       enddo


       write(*,*)
       write(*,*) 'WENO SCHEME(FINITE DIFFERENCE METHOD):'
       write(*,*) '****************************************'


       call exact(N,dx,X,t,u_exact)

       call WENO(N,dx,X,dt,nt,u_avg)

       !!! calculate error
       u_error = 0.0d0
       error_0 = 0.0d0
       do i=0,N-1
              u_error = ABS(u_exact(i)-u_avg(i))+u_error
              error_0 = MAX(error_0, ABS(u_exact(i)-u_avg(i)))
       enddo
       u_error = u_error/N

       write(*,*) 'Number of points is:',N
       write(*,*) 'alpha=', alpha
       write(*,*) 'L1 Normal is:', u_error
       !write(*,*) 'L Infinity Normal is:', error_0


       open(unit=01,file='./result.dat',status='unknown')
       write(01,101)
       write(01,102)
       write(01,103) N

       do i = 0,N-1
                     write(01,100) X(i), u_avg(i),u_exact(i)
       enddo


       close(01)
       print*,'****************************************'


100    format(2x,10(e12.6,'      '))
101    format('Title="Burgers Equation(WENO Scheme)"')
102    format('Variables=x,u_num,u_exact')
103    format('zone',1x,'i=',1x,i5,2x,'f=point')
       stop
       end program main



       subroutine exact(N,dx,X,t,u_exact)
       implicit none
       integer :: i, N, j, k
       real(8) :: dx, x0, temp, t
       real(8) :: X(0:N), u_exact(0:N)

       do i=0,N
              x0 = X(i)
              do
                     temp = x0
                     x0 = x0-(x0+t*(1.0d0/3.0d0+2.0d0/3.0d0*SIN(x0))-X(i))/(1.0d0+t*2.0d0/3.0d0*COS(x0))
                     if(ABS(temp-x0).GT.1d-15) then
                            cycle
                     else
                            exit
                     endif
              enddo
              u_exact(i) = 1.0d0/3.0d0+2.0d0/3.0d0*SIN(x0)
              !write(*,*) 'x0=',X(i),'u_exact=',u_exact(i)
       enddo



       return
       end subroutine exact



       subroutine WENO(N,dx,X,dt,nt,u_avg)
       implicit none
       integer :: N, i, j, nt
       real(8) :: dx, dt
       real(8) :: X(0:N), u_avg(0:N-1), du_avg(0:N-1), u1(0:N-1), u2(0:N-1)

       do i=0,N-1
              u_avg(i) = 1.0d0/3.0d0+2.0d0/3.0d0*SIN(X(i))
       enddo

       do j=1,nt
              call Lax(u_avg, du_avg, dx, N)
              do i=0,N-1
                     u1(i) = u_avg(i)+dt*du_avg(i)
              enddo

              call Lax(u1,du_avg,dx,N)
              do i=0,N-1
                     u2(i) = 3.0d0/4.0d0*u_avg(i)+1.0d0/4.0d0*(u1(i)+dt*du_avg(i))
              enddo

              call Lax(u2,du_avg,dx,N)
              do i=0,N-1
                     u_avg(i) = 1.0d0/3.0d0*u_avg(i)+2.0d0/3.0d0*(u2(i)+dt*du_avg(i))
              enddo
       enddo

!       do i=0,N-1
!              write(*,*) 'x0=',X(i),'u(i)=',u_avg(i)
!       enddo
!       write(*,*) 'x0=',X(N),'u(i)=',u_avg(0)

       return
       end subroutine WENO

       subroutine Recon(u_avg, u_reconl, u_reconr, dx, N)
       implicit none
       integer :: i, N, j
       real(8) :: u_avg(0:N-1), u_reconl(0:N-1), u_reconr(0:N-1), dx
       real(8) :: Beta(0:2), ul(0:2), ur(0:2)
       real(8) :: epsilon

       epsilon = 0.D0
       DO i=0, N-1
              epsilon = epsilon+ABS(u_avg(i))
       enddo
       epsilon = epsilon*MIN(1.D-8, dx*dx)/DBLE(N)
       do i=0,N-1

              ul(0) = 1.0d0/3.0d0*u_avg(MOD(i-2+N,N))-7.0d0/6.0d0*u_avg(MOD(i-1+N,N))+11.0d0/6.0d0*u_avg(i)
              ur(0) = -1.0d0/6.0d0*u_avg(MOD(i-2+N,N))+5.0d0/6.0d0*u_avg(MOD(i-1+N,N))+1.0d0/3.0d0*u_avg(i)

              ul(1) = -1.0d0/6.0d0*u_avg(MOD(i-1+N,N))+5.0d0/6.0d0*u_avg(i)+1.0d0/3.0d0*u_avg(MOD(i+1,N))
              ur(1) = -1.0d0/6.0d0*u_avg(MOD(i+1,N))+5.0d0/6.0d0*u_avg(i)+1.0d0/3.0d0*u_avg(MOD(i-1+N,N))

              ul(2) = 1.0d0/3.0d0*u_avg(i)+5.0d0/6.0d0*u_avg(MOD(i+1,N))-1.0d0/6.0d0*u_avg(MOD(i+2,N))
              ur(2) = 1.0d0/3.0d0*u_avg(MOD(i+2,N))-7.0d0/6.0d0*u_avg(MOD(i+1,N))+11.0d0/6.0d0*u_avg(i)

              Beta(0) = 13.d0/12.d0*(u_avg(MOD(i-2+N, N))-2.d0*u_avg(MOD(i-1+N, N))+u_avg(i))**2 &
                     +0.25d0*(u_avg(MOD(i-2+N, N))-4.0d0*u_avg(MOD(i-1+N, N))+3.d0*u_avg(i))**2
              Beta(1) = 13.d0/12.d0*(u_avg(MOD(i-1+N, N))-2.d0*u_avg(i)+u_avg(MOD(i+1, N)))**2+ &
                     0.25d0*(u_avg(MOD(i-1+N, N))-u_avg(MOD(i+1, N)))**2
              Beta(2) = 13.d0/12.d0*(u_avg(i)-2.d0*u_avg(MOD(i+1, N))+u_avg(MOD(i+2, N)))**2+ &
                     0.25d0*(3.d0*u_avg(i)-4.d0*u_avg(MOD(i+1, N))+u_avg(MOD(i+2, N)))**2
              Beta(0) = 0.1d0/((Beta(0)+epsilon)*(Beta(0)+epsilon))
              Beta(1) = 0.6d0/((Beta(1)+epsilon)*(Beta(1)+epsilon))
              Beta(2) = 0.3d0/((Beta(2)+epsilon)*(Beta(2)+epsilon))

              epsilon = 0.d0
              DO j=0, 2
                     epsilon = epsilon+Beta(j)
              enddo
              DO j=0, 2
                     Beta(j) = Beta(j)/epsilon
              enddo
              u_reconl(i) = Beta(0)*ul(0)+Beta(1)*ul(1)+Beta(2)*ul(2)
              u_reconr(MOD(i-1+N, N)) = Beta(0)*ur(2)+Beta(1)*ur(1)+Beta(2)*ur(0)
       enddo
       end subroutine Recon

       subroutine Lax(u_avg, du_avg, dx, N)
       implicit none
       integer :: i, N, s
       real(8) :: u_avg(0:N-1), du_avg(0:N-1), dx
       real(8) :: flux_P(0:N-1), flux_N(0:N-1)
       real(8) :: u_reconl(0:N-1), u_reconr(0:N-1), f(0:N-1)
       real(8) :: alpha

       alpha=0.d0
       do i=0, N-1
              alpha = MAX(alpha, ABS(u_avg(i)))
       enddo
       do i=0, N-1
              flux_P(i) = 0.5d0*(0.5d0*u_avg(i)*u_avg(i)+alpha*u_avg(i))
              flux_N(i) = 0.5*(0.5d0*u_avg(i)*u_avg(i)-alpha*u_avg(i))
       enddo
       call Recon(flux_P, u_reconl, f, dx, N)
       call Recon(flux_N, f, u_reconr, dx, N)
       do i=0, N-1
              f(i) = u_reconl(i)+u_reconr(i)
       enddo
       do i=0, N-1
              du_avg(i) = -(f(i)-f(MOD(i-1+N, N)))/dx
       enddo

       return
       end subroutine Lax
