
!!!    This program solves an ordinary diffential equation by using Euler Predictor-Corrector method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       real :: x0,xn,y0,h

       write(*,*) 'Euler Method: '
       write(*,*) 'ODE: y`=y^2-x^2'
       x0 = 0.0
       xn = 1.0
       y0 = 1.0
       h = 0.01
       write(*,*) 'B.C.: x(0)=',x0,'x(N)=',xn,'y(0)=',y0,'h=',h
       call Euler(x0,xn,y0,h)

       stop
       end program main



       subroutine Euler(x0,xn,y0,h)
       implicit none
       integer :: i,n
       real :: h,x0,xn,y0,yp,yc,y1,x1,f

       write(*,*)'x=',x0,'  y=',y0
       n=(xn-x0)/h
       i=1
       do
              x1=x0+h
              yp=y0+h*f(x0,y0)
              yc=y0+h*f(x1,yp)
              y1=(yp+yc)/2
              write(*,*)'x=',x1,'  y=',y1
              if(i.LT.n)then
                     x0=x1
                     y0=y1
                     i=i+1
              else
                     exit
              endif
       enddo
       end subroutine Euler



       function f(x,y)
       f = y*y - x*x
       end function f



