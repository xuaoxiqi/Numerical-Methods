

!!!    This program solves Burgers Equation using WENO Scheme in Finite Volume Formulation.
!!!    with Lax-Friedrichs Flux
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

        program main
        implicit none
        integer, parameter :: N=320
        real(8), parameter :: Pi=3.1415926535898d0
        integer :: i, nt
        real(8) :: dx, dt, t, error_0, error_1
        real(8) :: X(N), u_avg(-1:N+3), u_exact(N)

        dx = 2.0d0*Pi/float(N-1)
        dt = 0.01d0*dx
        t = 0.1d0*Pi
        nt = NINT(t/dt)
        do i=1,N
            X(i) = (i-1)*dx
        enddo

        write(*,*)
        write(*,*) 'WENO Scheme in Finite Volume Formulation:'
        write(*,*) 'with Lax-Friedrichs Flux'
        write(*,*) '*********************************************'

        call exact(N,dx,X,t,u_exact)

        call WENO(N,dx,X,dt,nt,u_avg)

       !!! calculate error
        error_0 = 0.0d0
        error_1 = 0.0d0
        do i=1,N-1
            error_0 = MAX(error_0, ABS(u_exact(i)-u_avg(i)))
            error_1 = ABS(u_exact(i)-u_avg(i))+error_1
        enddo
        error_1 = error_1/float(N-1)

        write(*,*) 'N =',N
        write(*,*) 't =',t
        write(*,*) 'dx =',dx
        write(*,*) 'dt =',dt
        write(*,*) 'LInfinity Normal is:', error_0
        write(*,*) 'L1        Normal is:', error_1

        open(unit=01,file='./result.dat',status='unknown')
        write(01,101)
        write(01,102)
        write(01,103) N-1

        do i = 1,N-1
            write(01,100) X(i)+0.5d0*dx, u_avg(i),u_exact(i)
        enddo

        close(01)

100     format(2x,10(e12.6,'      '))
101     format('Title="Burgers Equation(WENO Scheme)"')
102     format('Variables=x,u_num,u_exact')
103     format('zone',1x,'i=',1x,i5,2x,'f=point')

        stop
        end program main



        subroutine exact(N,dx,X,t,u_exact)
        implicit none
        integer :: i, N, j, k
        real(8) :: dx, x0, temp, t
        real(8) :: X(N)
        real(8) :: u_exact(N), u_t(10), u_x(10), weight(10)
        data u_t/0.9931285991850949d0,0.9639719272779138d0,0.9122344282513259d0,0.8391169718222188d0,0.7463319064601508d0, &
              0.6360536807265150d0,0.5108670019508271d0,0.3737060887154196d0,0.2277858511416451d0,0.07652652113349734d0/  !!! x_k
        data weight/0.01761400713915212d0,0.04060142980038694d0,0.06267204833410906d0,0.08327674157670475d0,0.1019301198172404d0, &
              0.1181945319615184d0,0.1316886384491766d0,0.1420961093183821d0,0.1491729864726037d0,0.1527533871307259d0/  !!! weight coefficient A_k

        do i=1,N
            u_exact(i)=0.d0
            do k=-1,1,2

                do j=1,10
                    u_x(j) = X(i)+0.5d0*dx*float(k)*u_t(j)
                enddo

                do j=1,10
                    x0 = u_x(j)
                    do
                        temp = x0
                        x0 = x0-(x0+t*(1.0d0/3.0d0+2.0d0/3.0d0*SIN(x0))-u_x(j))/(1.0d0+t*2.0d0/3.0d0*COS(x0))
                        if(ABS(temp-x0).GT.1e-15) then
                            cycle
                        else
                            exit
                        endif
                    enddo
                    u_exact(i) = u_exact(i)+weight(j)*(1.0d0/3.0d0+2.0d0/3.0d0*SIN(x0)) *0.5d0
                enddo

            enddo
        enddo

        return
        end subroutine exact


        subroutine WENO(N,dx,X,dt,nt,u_avg)
        implicit none
        integer :: N, i, j, nt
        real(8) :: dx, dt
        real(8) :: X(N), u_avg(-1:N+3), u1(-1:N+3), u2(-1:N+3), du_avg(N)

        do i=1,N
            u_avg(i) = 1.0d0/3.0d0+4.0d0/3.0d0/dx*SIN(X(i))*SIN(0.5d0*dx)
        enddo
        u_avg(N+1) = u_avg(2)
        u_avg(N+2) = u_avg(3)
        u_avg(N+3) = u_avg(4)
        u_avg(0) = u_avg(N-1)
        u_avg(-1) = u_avg(N-2)

        do j=1,nt

            call Lax(dx,N,u_avg,du_avg)
            do i=1,N
                u1(i) = u_avg(i)+dt*du_avg(i)
            enddo
            u1(N+1) = u1(2)
            u1(N+2) = u1(3)
            u1(N+3) = u1(4)
            u1(0) = u1(N-1)
            u1(-1) = u1(N-2)

            call Lax(dx,N,u1,du_avg)
            do i=1,N
                u2(i) = 3.0d0/4.0d0*u_avg(i)+1.0d0/4.0d0*(u1(i)+dt*du_avg(i))
            enddo
            u2(N+1) = u2(2)
            u2(N+2) = u2(3)
            u2(N+3) = u2(4)
            u2(0) = u2(N-1)
            u2(-1) = u2(N-2)

            call Lax(dx,N,u2,du_avg)
            do i=1,N
                u_avg(i) = 1.0d0/3.0d0*u_avg(i)+2.0d0/3.0d0*(u2(i)+dt*du_avg(i))
            enddo
            u_avg(N+1) = u_avg(2)
            u_avg(N+2) = u_avg(3)
            u_avg(N+3) = u_avg(4)
            u_avg(0) = u_avg(N-1)
            u_avg(-1) = u_avg(N-2)

        enddo

        return
        end subroutine WENO


        subroutine Lax(dx,N,u_avg,du_avg)
        implicit none
        integer :: i, N
        real(8) :: u_avg(-1:N+3), du_avg(N), dx
        real(8) :: flux_P(0:N-1), flux_N(0:N-1)
        real(8) :: u_l(N+1), u_r(N+1), f(0:N)
        real(8) :: alpha

        alpha=0.0d0
        do i=1,N
            alpha = MAX(alpha, ABS(u_avg(i)))
        enddo

        call Recon(dx,N,u_avg,u_l,u_r)

        do i=1,N
            f(i) = 0.5d0*(0.5d0*u_l(i+1)*u_l(i+1)+0.5d0*u_r(i)*u_r(i)-alpha*(u_r(i)-u_l(i+1)))
        enddo
        f(0) = f(N-1)

        do i=1,N
            du_avg(i) = -(f(i)-f(i-1))/dx
        enddo

        return
        end subroutine Lax


        subroutine Recon(dx,N,u_avg,u_l,u_r)
        implicit none
        integer :: i, N, j
        real(8) :: u_avg(-1:N+3), u_l(N+1), u_r(N+1), dx
        real(8) :: beta(0:2)
        real(8) :: alphal(0:2), alphar(0:2), omgl(0:2), omgr(0:2), ul(0:2), ur(0:2)
        real(8) :: epsilon

        epsilon = 1e-6

        do i=1,N+1

            ur(0) = 1.0d0/3.0d0*u_avg(i)+5.0d0/6.0d0*u_avg(i+1)-1.0d0/6.0d0*u_avg(i+2)
            ur(1) = -1.0d0/6.0d0*u_avg(i-1)+5.0d0/6.0d0*u_avg(i)+1.0d0/3.0d0*u_avg(i+1)
            ur(2) = 1.0d0/3.0d0*u_avg(i-2)-7.0d0/6.0d0*u_avg(i-1)+11.0d0/6.0d0*u_avg(i)

            ul(0) = 11.0d0/6.0d0*u_avg(i)-7.0d0/6.0d0*u_avg(i+1)+1.0d0/3.0d0*u_avg(i+2)
            ul(1) = 1.0d0/3.0d0*u_avg(i-1)+5.0d0/6.0d0*u_avg(i)-1.0d0/6.0d0*u_avg(i+1)
            ul(2) = -1.0d0/6.0d0*u_avg(i-2)+5.0d0/6.0d0*u_avg(i-1)+1.0d0/3.0d0*u_avg(i)

            !compute weights
            beta(0) = 13.d0/12.d0*(u_avg(i)-2.0d0*u_avg(i+1)+u_avg(i+2))**2+ &
                     0.25d0*(3.d0*u_avg(i)-4.d0*u_avg(i+1)+u_avg(i+2))**2
            beta(1) = 13.d0/12.d0*(u_avg(i-1)-2.d0*u_avg(i)+u_avg(i+1))**2+ &
                     0.25d0*(u_avg(i-1)-u_avg(i+1))**2
            beta(2) = 13.0d0/12.0d0*(u_avg(i-2)-2.d0*u_avg(i-1)+u_avg(i))**2 &
                     +0.25d0*(u_avg(i-2)-4.0d0*u_avg(i-1)+3.d0*u_avg(i))**2

            alphar(0) = 0.3d0/(beta(0)+epsilon)/(beta(0)+epsilon)
            alphar(1) = 0.6d0/(beta(1)+epsilon)/(beta(1)+epsilon)
            alphar(2) = 0.1d0/(beta(2)+epsilon)/(beta(2)+epsilon)

            alphal(0) = 0.1d0/(beta(0)+epsilon)/(beta(0)+epsilon)
            alphal(1) = 0.6d0/(beta(1)+epsilon)/(beta(1)+epsilon)
            alphal(2) = 0.3d0/(beta(2)+epsilon)/(beta(2)+epsilon)

            do j=0, 2
                omgr(j) = alphar(j)/(alphar(0)+alphar(1)+alphar(2))
                omgl(j) = alphal(j)/(alphal(0)+alphal(1)+alphal(2))
            enddo

            u_r(i) = omgr(0)*ur(0)+omgr(1)*ur(1)+omgr(2)*ur(2)
            u_l(i) = omgl(0)*ul(0)+omgl(1)*ul(1)+omgl(2)*ul(2)

        enddo

        return
        end subroutine Recon

