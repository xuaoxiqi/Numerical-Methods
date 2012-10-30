
!!!    This program solves Burgers Equation using ENO Scheme.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       integer, parameter :: N=160
       real(8), parameter :: Pi=3.1415926535897932385d0
       integer :: i, nt
       real(8) :: alpha, dx, dt, t, u_error, error_0
       real(8) :: u_avg(0:N-1)
       real(8) :: u_exact(0:N), X(0:N)

       alpha = 0.5d0
       dx = 2.0d0*Pi/float(N)
       dt = alpha*dx
       t = 0.1d0*Pi
       nt = NINT(t/dt)

       do i=0,N
              X(i) = DBLE(i)*dx
       enddo


       write(*,*) ''
       write(*,*) 'ENO SCHEME:'
       write(*,*) '****************************************'


       call exact(N,dx,X,t,u_exact)

       call ENO(N,dx,X,dt,nt,u_avg)

       !!! calculate error
       u_error = 0.0d0
       error_0 = 0.0d0
       do i=0,N-1
              u_error = ABS(u_exact(i)-u_avg(i))+u_error
              error_0 = MAX(error_0, ABS(u_exact(i)-u_avg(i)))
       enddo
       u_error = u_error/DBLE(N)

       write(*,*) 'Number of points is:',N
       write(*,*) 'alpha =',alpha
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
101    format('Title="Burgers Equation(ENO Scheme)"')
102    format('Variables=x,u_num,u_exact')
103    format('zone',1x,'i=',1x,i5,2x,'f=point')
       stop
       end program main



       subroutine exact(N,dx,X,t,u_exact)
       implicit none
       integer :: i, N, j, k
       real(8) :: dx, x0, temp, t
       real(8) :: X(0:N)
       real(8) :: u_exact(0:N), u_t(10), u_x(10), weight(10)
       data u_t/0.9931285991850949d0,0.9639719272779138d0,0.9122344282513259d0,0.8391169718222188d0,0.7463319064601508d0, &
              0.6360536807265150d0,0.5108670019508271d0,0.3737060887154196d0,0.2277858511416451d0,0.07652652113349734d0/  !!! x_k
       data weight/0.01761400713915212d0,0.04060142980038694d0,0.06267204833410906d0,0.08327674157670475d0,0.1019301198172404d0, &
              0.1181945319615184d0,0.1316886384491766d0,0.1420961093183821d0,0.1491729864726037d0,0.1527533871307259d0/  !!! weight coefficient A_k

       do i=0,N
              u_exact(i)=0.d0
              do k=-1,1,2
                     do j=1,10
                     u_x(j) = X(i)-dx*0.5d0+dx*(1.0d0+DBLE(k)*u_t(j))/2.d0
                     enddo
                     do j=1,10
                            x0 = u_x(j)
                            do
                                   temp = x0
                                   x0 = x0-(x0+t*(1.0d0/3.0d0+2.0d0/3.0d0*SIN(x0))-u_x(j))/(1.0d0+t*2.0d0/3.0d0*COS(x0))
                                   if(ABS(temp-x0).GT.1d-14) then
                                          cycle
                                   else
                                          exit
                                   endif
                            enddo
                            u_exact(i) = u_exact(i)+weight(j)*(1.0d0/3.0d0+2.0d0/3.0d0*SIN(x0))*0.5d0
                     enddo
              enddo

              !write(*,*) 'x0=',X(i),'u_exact=',u_exact(i)
       enddo



       return
       end subroutine exact



        subroutine ENO(N,dx,X,dt,nt,u_avg)
        implicit none
        integer :: N, i, j, nt
        real(8) :: dx, dt, ul, ur
        real(8) :: X(0:N), u_avg(0:N-1), du_avg(0:N-1), f(N), u1(0:N-1), u2(0:N-1)

        do i=1,N-1
            u_avg(i) = 1.0d0/3.0d0+4.0d0/3.0d0/dx*SIN(X(i))*SIN(0.5d0*dx)
        enddo
        u_avg(0) = 1.0d0/3.0d0

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

        return
        end subroutine ENO


        subroutine Godunov(u_avg, du_avg, dx, N)
        implicit none
        integer :: i, N, s
        real(8) :: u_avg(0:N-1), du_avg(0:N-1), dx
        real(8) :: u_reconl(0:N-1), u_reconr(0:N-1), f(0:N-1)

        call Recon(u_avg,u_reconl,u_reconr,dx,N)
        do i=0, N-1
            if(u_reconl(i).LE.u_reconr(i)) then
                f(i) = MIN(0.5d0*u_reconl(i)*u_reconl(i),0.5d0*u_reconr(i)*u_reconr(i))
                if(u_reconl(i)*u_reconr(i).LT.0.0d0) f(i) = 0.0d0
            else
                f(i) = MAX(0.5d0*u_reconl(i)*u_reconl(i),0.5d0*u_reconr(i)*u_reconr(i))
            endif
        enddo

        do i=0,N-1
            du_avg(i) = -(f(i)-f(MOD(i-1+N, N)))/dx
        enddo

        return
        endsubroutine Godunov


        subroutine Recon(u_avg, u_reconl, u_reconr, dx, N)
        implicit none
        integer :: i, N, s
        real(8) :: u_avg(0:N-1), u_reconl(0:N-1), u_reconr(0:N-1), dx

        do i=0,N-1
            call stencil(u_avg,3,i,N,s)
            select case (s)
                case (-2)
                    u_reconl(i) = 1.0d0/3.0d0*u_avg(MOD(i-2+N,N))-7.0d0/6.0d0*u_avg(MOD(i-1+N,N))+11.0d0/6.0d0*u_avg(i)
                    u_reconr(MOD(i-1+N, N)) = -1.0d0/6.0d0*u_avg(MOD(i-2+N,N))+5.0d0/6.0d0*u_avg(MOD(i-1+N,N))+1.0d0/3.0d0*u_avg(i)
                case (0)
                    u_reconl(i) = -1.0d0/6.0d0*u_avg(MOD(i-1+N,N))+5.0d0/6.0d0*u_avg(i)+1.0d0/3.0d0*u_avg(MOD(i+1,N))
                    u_reconr(MOD(i-1+N, N)) = -1.0d0/6.0d0*u_avg(MOD(i+1,N))+5.0d0/6.0d0*u_avg(i)+1.0d0/3.0d0*u_avg(MOD(i-1+N,N))
                case (2)
                    u_reconl(i) = 1.0d0/3.0d0*u_avg(i)+5.0d0/6.0d0*u_avg(MOD(i+1,N))-1.0d0/6.0d0*u_avg(MOD(i+2,N))
                    u_reconr(MOD(i-1+N, N)) = 1.0d0/3.0d0*u_avg(MOD(i+2,N))-7.0d0/6.0d0*u_avg(MOD(i+1,N))+11.0d0/6.0d0*u_avg(i)
            end select
        enddo

        return
        end subroutine Recon


       subroutine stencil(u, m, k, N, s)
        implicit none
        integer :: i, j, m, s, k, N, k_temp
        real(8) :: u(0:N-1), u_differ(0:N-1)

        do j=0, N-1
            u_differ(j) = u(MOD(j+1,N))-u(j)
        enddo

        k_temp = k;
        s = 0

        do i=2, m

            if(u_differ(MOD(k_temp-1+N,N)).LE.u_differ(MOD(k_temp+N,N))) then
                s = s-1
                k_temp = k_temp-1
            else
                s = s+1
            endif

            do j=0, N-1
                u_differ(j) = u_differ(MOD(j+1,N))-u_differ(j)
            enddo

        enddo

        return
        end subroutine stencil

