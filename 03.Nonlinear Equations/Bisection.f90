
!!!    This program finds all the roots of nonlinear equation by using Bisection Method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

    program main
    implicit none
    integer, parameter ::dp=kind(0.d0)
    real(dp) :: a0,b0,h

    write(*,*) 'Bisection Method:'
    write(*,*) 'Eq: 2*exp(-x)-sin(x)=0'
    a0 = 0.0_dp
    b0 = 10.0_dp
    h = 0.02_dp
    write(*,*) 'Roots in interval[',a0,b0,'] step h=',h
    call Bisection(a0,b0,h)

    stop
    end program main



    subroutine Bisection(a0,b0,h)
    implicit none
    integer, parameter ::dp=kind(0.d0)
    real(dp) :: a0,b0,h,z,a,b,c,f

    z = a0
    do while(z.LT.b0)
        if(ABS(f(z)).LT.0.01_dp) then
            write(*,*) 'x=',z
        else
            if(f(z)*f(z+h).LT.0.0_dp) then
                a = z
                b = z+h
                c = (a+b)/2.0_dp
                if(ABS(f(c)).LT.0.01_dp) then
                    write(*,*) 'x=',c
                else
                    if(f(a)*f(c).LT.0.0_dp) b = c
                    if(f(b)*f(c).LT.0.0_dp) a = c
                endif
                else
                    continue
                endif
            endif
        z=z+h
    enddo

    return
    end subroutine Bisection



    function f(x)
    implicit none
    integer, parameter ::dp=kind(0.d0)
    real(dp) :: f, x
    f = 2.0_dp*exp(-x)-sin(x)
    end function f
