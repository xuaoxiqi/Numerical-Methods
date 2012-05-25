
!!!    This program integrates function f(x) by using Simpson rule with doubling number of intervals.
!!!    The result is modified by Cotes Method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       real(8) :: a, b, f, integral, eps
       external f

       a = 0.0
       b = 2.0
       eps = 1.0e-6
       write(*,*) 'Simpson Method(result modified by Cotes Method):'
       write(*,*) 'f(x) = sin(0.5*Pi*x^2)'
       write(*,*) 'Integration on [',a,b,']'
       call Simpson(f,a,b,eps,integral)
       write(*,*) 'Result:',integral

       stop
       end program main



       subroutine Simpson(f,a,b,eps,integral)
       implicit none
       real(8) :: f, a, b, eps, integral
       real(8) :: sn, s2n, h, x
       real(8), parameter :: coeff = 1.0/15.0 ! error estimate coeff
       integer, parameter :: nmax=1048576              ! max number of intervals
       integer n, i

       ! evaluate integral for 2 intervals (three points)
       h = (b-a)/2.0
       sn = (1.0/3.0)*h*(f(a)+4.0*f(a+h)+f(b))
       ! loop over number of intervals (starting from 4 intervals)
       n=4
       do while (n <= nmax)
              s2n = 0.0
              h = (b-a)/dfloat(n)
              do i=2, n-2, 2
                     x = a+dfloat(i)*h
                     s2n = s2n + 2.0*f(x) + 4.0*f(x+h)
              enddo
              s2n = (s2n + f(a) + f(b) + 4.0*f(a+h))*h/3.0
!!!           Similar as Cotes Method
              if(coeff*ABS(s2n-sn) <= eps) then
                     integral = s2n + coeff*(s2n-sn)
                     exit
              end if
              sn = s2n
              n = n*2
       enddo

       return
       end subroutine Simpson



       function f(x)
       implicit none
       real(8), parameter :: Pi=3.141592653589793
       real(8) :: f, x

       f=SIN(Pi*x*x/2.0)

       return
       end function f


