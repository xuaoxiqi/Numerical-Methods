
!!!    This program solves an ordinary diffential equation by using 4-Order Runge-Kutta method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

    program main
    implicit none
    integer, parameter :: N=101
    integer, parameter :: dp=kind(0.d0)
    real(dp) :: x(N), y(N)
    real(dp) :: dx
    integer :: i

    dx = (1.0_dp-0.0_dp)/(N-1)
    do i=1,N
        x(i) = (i-1)*dx
    enddo
    y = 0.0_dp
    y(1) = 1.0_dp

    write(*,*) 'Euler Method: '
    write(*,*) 'ODE: y`=y^2-x^2'
    write(*,*) 'dx = ',dx
    write(*,*) 'B.C.: x(1)=',x(1),'x(N)=',x(N)
    write(*,*) 'B.C.: y(1)=',y(1)

    call RK4(x,y,dx,N)

    stop
    end program main



    subroutine RK4(x,y,dx,N)
    implicit none
    integer :: i, N
    integer, parameter :: dp=kind(0.d0)
    real(dp) :: dx, x(N), y(N)
    real(dp) :: k1, k2, k3, k4, f

    do i=1,N-1
        k1 = f(x(i),y(i))
        k2 = f((x(i)+0.5_dp*dx),(y(i)+0.5_dp*dx*k1))
        k3 = f((x(i)+0.5_dp*dx),(y(i)+0.5_dp*dx*k2))
        k4 = f(x(i+1),(y(i)+dx*k3))
        y(i+1) = y(i)+dx*(k1+2.0_dp*k2+2.0_dp*k3+k4)/6.0_dp
        write(*,*) 'x=',x(i+1),'y=',y(i+1)
    enddo

    return
    end subroutine RK4



    function f(x,y)
    integer, parameter:: dp=kind(0.d0)
    real(dp) :: f, x, y
    f = y*y - x*x
    end function f
