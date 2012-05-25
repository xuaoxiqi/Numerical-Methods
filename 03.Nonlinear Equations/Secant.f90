
!!!    This program finds all the roots of nonlinear equation by using Secant Method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       real :: a, b, root, eps

       write(*,*) 'Secant Method:'
       write(*,*) 'Eq: 2*exp(-x)-sin(x)=0'
       a=2.0
       b=4.0
       eps = 1.0e-5
       write(*,*) 'Root closely near interval[',a,b,']'
       call Secant(a,b,root,eps)
       write(*,*) 'x=',root

       stop
       end program main



       subroutine Secant(x1,x2,root,eps)
       implicit none
       integer, parameter:: itc=200
       integer :: i
       real :: x1, x2, x3, df, root, eps,f

!!!    Iterative refining the solution
       do i=1,itc
              df = (x2-x1)/(f(x2)-f(x1))
              x3 = x2 - f(x2)*df
!!!           check the step. if it is improbably large - use bisection
              if(abs(x3) > 100.0*abs(x2)) x3 = (x2+x1)/2.0
!!!           condition to stop iterations
              if(abs(f(x3))<= eps) exit
              x1 = x2;
              x2 = x3;
       end do
       root=x3

       return
       end subroutine Secant



       function f(x)
       f=2*exp(-x)-sin(x)
       end
