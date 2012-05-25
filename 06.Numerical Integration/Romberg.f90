
!!!    This program integrates function f(x) by using Romberg Method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       real(8) :: f, a, b, eps, x, ri
       external f

       a = 0.0
       b = 5.0
       eps = 1.0e-6
       write(*,*) 'Romberg Method:'
       write(*,*) 'f(x) = cos(0.5*Pi*x*x)'
       write(*,*) 'Integration on [',a,b,']'
       call Romberg(a,b,eps,x,ri)
       write(*,*) 'Result:',ri

       stop
       end program main



       subroutine Romberg(a,b,eps,x,ri)
       implicit none
       integer :: k
       real(8) ::  a, b, eps, x, ri, ri0, f, h, fa, fb, t1, t2, s, s1, s2, c1, c2, ff

       h = b-a
       fa = f(a)
       fb = f(b)
       t1 = h*(fa+fb)/2.0
       k = 1
       do
              s = 0.0
              x = a+0.5*h
              do while(x.LT.b)
                     ff = f(x)
                     s = s+ff
                     x = x+h
              enddo
              t2 = (t1+h*s)/2.0
              s2 = t2+(t2-t1)/3.0
              if(k.NE.1) then
                     c2 = s2+(s2-s1)/15.0
                     if(k.NE.2) then
                            ri = c2+(c2-c1)/63.0
                            if(k.NE.3) then
                                   if(ABS(ri-ri0).GT.ABS(ri)*eps) then
                                          ri0 = ri
                                          c1 = c2
                                          k = k+1
                                          h = 0.5*h
                                          t1 = t2
                                          s1 = s2
                                          cycle
                                   else
                                          exit
                                   endif
                            else
                                   ri0 = ri
                                   c1 = c2
                                   k = k+1
                                   h = 0.5*h
                                   t1 = t2
                                   s1 = s2
                                   cycle
                            endif

                     else
                            c1 = c2
                            k = k+1
                            h = 0.5*h
                            t1 = t2
                            s1 = s2
                            cycle
                     endif

              else
                     k = k+1
                     h = h/2.0
                     t1 = t2
                     s1 = s2
                     cycle
              endif
       enddo


       return
       end subroutine Romberg



       function f(x)
       implicit none
       real(8), parameter :: Pi=3.141592653589793
       real(8) :: f, x

       f=COS(0.5*Pi*x*x) !!! Function for integration

       return
       end function f
