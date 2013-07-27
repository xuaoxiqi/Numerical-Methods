
!!!    This program solves an ordinary diffential equation by using Euler Predictor-Corrector method.
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

    call Euler(x,y,dx,N)

    stop
    end program main



    subroutine Euler(x,y,dx,N)
    implicit none
    integer, parameter :: dp=kind(0.d0)
    real(dp) :: dx, x(N), y(N)
    real(dp) :: yp, yc, f
    integer :: i, N

    do i=1,N-1
        yp = y(i)+dx*f(x(i),y(i))
        yc = y(i)+dx*f(x(i+1),yp)
        y(i+1) = (yp+yc)/2.0_dp
        write(*,*) "x=",x(i+1), "y=", y(i+1)
    enddo

    return
    end subroutine Euler



    function f(x,y)
    implicit none
    integer, parameter:: dp=kind(0.d0)
    real(dp) :: f, x, y
    f = y*y - x*x
    end function f

