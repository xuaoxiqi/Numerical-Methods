

!!!    This program solves Burgers Equation using WENO Scheme with Finite Difference Method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

        program main
        implicit none
        integer, parameter :: N=10
        real(8), parameter :: Pi=3.1415926535897932385d0
        integer :: i, nt
        real(8) :: dx, dt, t, error_0, error_1
        real(8) :: alpha
        real(8) :: u(0:N-1)
        real(8) :: u_exact(0:N), X(0:N-1)

        dx = 2.0d0*Pi/float(N)
        !!!dt = 0.1*dx**(5.0d0/3.0d0)
        alpha = 0.01d0
        dt = alpha*dx
        t = 0.1*Pi
        nt = NINT(t/dt)

        do i=0,N-1
            X(i) = float(i)*dx
        enddo

        write(*,*)
        write(*,*) 'WENO SCHEME(FINITE DIFFERENCE METHOD):'
        write(*,*) '****************************************'

        call exact(N,dx,X,t,u_exact)

        call WENO(N,dx,X,dt,nt,u)

       !!! calculate error
        error_0 = 0.0d0
        error_1 = 0.0d0
        do i=0,N-1
            error_0 = MAX(error_0, ABS(u_exact(i)-u(i)))
            error_1 = ABS(u_exact(i)-u(i))+error_1
        enddo
        error_1 = error_1/float(N)

        write(*,*) 'Number of points is:',N
        write(*,*) 'dx =',dx
        write(*,*) 'dt =',dt
        !write(*,*) 'LInfinity Normal is:', error_0
        write(*,*) 'L1        Normal is:', error_1

        open(unit=01,file='./result.dat',status='unknown')
        write(01,101)
        write(01,102)
        write(01,103) N

        do i = 0,N-1
            write(01,100) X(i), u(i), u_exact(i)
        enddo

        close(01)
        print*,'****************************************'

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
        real(8) :: X(0:N-1), u_exact(0:N)

        do i=0,N-1
            x0 = X(i)
            do
                temp = x0
                x0 = x0-(x0+t*(1.0d0/3.0d0+2.0d0/3.0d0*SIN(x0))-X(i))/(1.0d0+t*2.0d0/3.0d0*COS(x0))
                if(ABS(temp-x0).GT.1e-15) then
                    cycle
                else
                    exit
                endif
            enddo
            u_exact(i) = 1.0d0/3.0d0+2.0d0/3.0d0*SIN(x0)
        enddo

        return
        end subroutine exact


        subroutine WENO(N,dx,X,dt,nt,u)
        implicit none
        integer :: N, i, j, nt
        real(8) :: dx, dt
        real(8) :: X(0:N-1), u(0:N-1), du(0:N-1), u1(0:N-1), u2(0:N-1)

        do i=0,N-1
            u(i) = 1.0d0/3.0d0+2.0d0/3.0d0*SIN(X(i))
        enddo

        do j=1,nt

            call Lax(dx,N,u,du)
            do i=0,N-1
                u1(i) = u(i)+dt*du(i)
            enddo

            call Lax(dx,N,u1,du)
            do i=0,N-1
                u2(i) = 3.0d0/4.0d0*u(i)+1.0d0/4.0d0*(u1(i)+dt*du(i))
            enddo

            call Lax(dx,N,u2,du)
            do i=0,N-1
                u(i) = 1.0d0/3.0d0*u(i)+2.0d0/3.0d0*(u2(i)+dt*du(i))
            enddo

        enddo


        return
        end subroutine WENO


        subroutine Lax(dx,N,u,du)
        implicit none
        integer :: i, N
        real(8) :: u(0:N-1), du(0:N-1), dx
        real(8) :: f_positive(0:N-1), f_negative(0:N-1)
        real(8) :: f_l(0:N-1), f_r(0:N-1), f(0:N-1)
        real(8) :: alpha

        alpha=0.0d0
        do i=0,N-1
            alpha = MAX(alpha, ABS(u(i)))
        enddo

        do i=0,N-1
            f_positive(i) = 0.5d0*(0.5d0*u(i)*u(i)+alpha*u(i))
        enddo
        call Recon(dx,N,f_positive,f,f_r)

        do i=0,N-1
            f_negative(i) = 0.5d0*(0.5d0*u(i)*u(i)-alpha*u(i))
        enddo
        call Recon(dx,N,f_negative,f_l,f)

        do i=0,N-1
            f(i) = f_l(MOD(i+1,N))+f_r(i)
        enddo

        do i=0,N-1
            du(i) = -(f(i)-f(MOD(i-1+N,N)))/dx
        enddo

        return
        end subroutine Lax


        !u_r = u_{i+1/2}
        !u_l = u_{i-1/2}

        subroutine Recon(dx,N,u,u_l,u_r)
        implicit none
        integer :: i, N, j
        real(8) :: u(0:N-1), u_l(0:N-1), u_r(0:N-1), dx
        real(8) :: beta(0:2)
        real(8) :: alphal(0:2), alphar(0:2), omgl(0:2), omgr(0:2), ul(0:2), ur(0:2)
        real(8) :: epsilon

        epsilon = 1e-6

        do i=0,N-1

            ur(0) = 1.0d0/3.0d0*u(i)+5.0d0/6.0d0*u(MOD(i+1,N))-1.0d0/6.0d0*u(MOD(i+2,N))
            ur(1) = -1.0d0/6.0d0*u(MOD(i-1+N,N))+5.0d0/6.0d0*u(i)+1.0d0/3.0d0*u(MOD(i+1,N))
            ur(2) = 1.0d0/3.0d0*u(MOD(i-2+N,N))-7.0d0/6.0d0*u(MOD(i-1+N,N))+11.0d0/6.0d0*u(i)

            ul(0) = 11.0d0/6.0d0*u(i)-7.0d0/6.0d0*u(MOD(i+1,N))+1.0d0/3.0d0*u(MOD(i+2,N))
            ul(1) = 1.0d0/3.0d0*u(MOD(i-1+N,N))+5.0d0/6.0d0*u(i)-1.0d0/6.0d0*u(MOD(i+1,N))
            ul(2) = -1.0d0/6.0d0*u(MOD(i-2+N,N))+5.0d0/6.0d0*u(MOD(i-1+N,N))+1.0d0/3.0d0*u(i)

            !compute weights
            beta(0) = 13.d0/12.d0*(u(i)-2.0d0*u(MOD(i+1,N))+u(MOD(i+2,N)))**2 &
                    +0.25d0*(3.d0*u(i)-4.d0*u(MOD(i+1,N))+u(MOD(i+2,N)))**2
            beta(1) = 13.d0/12.d0*(u(MOD(i-1+N,N))-2.0d0*u(i)+u(MOD(i+1,N)))**2 &
                    +0.25d0*(u(MOD(i-1+N,N))-u(MOD(i+1,N)))**2
            beta(2) = 13.0d0/12.0d0*(u(MOD(i-2+N,N))-2.d0*u(MOD(i-1+N,N))+u(i))**2 &
                    +0.25d0*(u(MOD(i-2+N,N))-4.0d0*u(MOD(i-1+N,N))+3.0d0*u(i))**2

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
