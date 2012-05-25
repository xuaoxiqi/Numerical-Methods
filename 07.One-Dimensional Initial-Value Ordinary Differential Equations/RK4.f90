
!!!    This program solves an ordinary diffential equation by using 4-Order Runge-Kutta method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       real :: x0,xn,y0,h

       write(*,*) '4-Order Runge-Kutta Method: '
       write(*,*) 'ODE: y`=y^2-x^2'
       x0 = 0.0
       xn = 1.0
       y0 = 1.0
       h = 0.1
       write(*,*) 'B.C.: x(0)=',x0,'x(N)=',xn,'y(0)=',y0,'h=',h

       call RK4(x0,xn,y0,h)

       stop
       end program main



       subroutine RK4(x0,xn,y0,h)
       implicit none
       integer :: i,n
       real :: h,x0,xn,y0,y1,x1,k1,k2,k3,k4,f

       write(*,*)'x=',x0,'   y=',y0
       n=(xn-x0)/h
       i=1
       do
              x1=x0+h
!!!    Get four sample values of the derivative.
              k1=f(x0,y0)
              k2=f((x0+0.5*h),(y0+0.5*h*k1))
              k3=f((x0+0.5*h),(y0+0.5*h*k2))
              k4=f(x1,(y0+h*k3))
!!!    Combine them to estimate the solution y1 at x1=x0+h
              y1=y0+h*(k1+2*k2+2*k3+k4)/6.0
              write(*,*)'x=',x1,'   y=',y1
              x0=x1
              y0=y1
              if(i.LT.n) then
                     i=i+1
              else
                     exit
              endif
       enddo

       return
       end subroutine RK4



       function f(x,y)
       f = y*y - x*x
       end function f





