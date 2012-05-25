
!!!    This program finds all the roots of nonlinear equation by using Newton's interation Method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       real :: x0

       write(*,*) 'Newton`s interation Method:'
       write(*,*) 'Eq: 2*exp(-x)-sin(x)=0'
       x0=6.0
       write(*,*) 'Root close near point x=',x0
       call Newton(x0)

       stop
       end program main



       subroutine Newton(x0)
       implicit none
       real :: x,x0,f,df

       x=x0
       x0=x-f(x)/df(x)
       do while(ABS(1-x/x0).GT.0.001)
              x=x0
              x0=x-f(x)/df(x)
       enddo
       write(*,*) 'x=',x0

       return
       end subroutine Newton



       function f(x)
       f=2*exp(-x)-sin(x)
       end function f

       function df(x)
       df=-2*exp(-x)-cos(x)
       end function df
