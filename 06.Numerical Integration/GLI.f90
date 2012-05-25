
!!!    This program integrates function f(x) with respect to x from a to b
!!!    by using Gauss-Legendre numerical integration method(n=20)
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       real(8) :: f, gl, a, b
       external f

       a = 0.0
       b = 1.0
       write(*,*) 'Gauss-Legendre Method(N=2*10):'
       write(*,*) 'f(x) = 1.0/(1.0+x)'
       write(*,*) 'Integration on [',a,b,']'
       call GLI(a,b,gl)
       write(*,*) 'Result:',gl

       stop
       end program main



       subroutine GLI(a,b,gl)
       implicit none
       integer, parameter :: n=10
       integer :: i
       real(8) :: t(n),w(n)
       real(8) :: a, b, c0, c1, c2, t1, t2, f, f1, f2, gl
       data t/0.9931285991850949,0.9639719272779138,0.9122344282513259,0.8391169718222188,0.7463319064601508, &
              0.6360536807265150,0.5108670019508271,0.3737060887154196,0.2277858511416451,0.07652652113349734/  !!! x_k
       data w/0.01761400713915212,0.04060142980038694,0.06267204833410906,0.08327674157670475,0.1019301198172404, &
              0.1181945319615184,0.1316886384491766,0.1420961093183821,0.1491729864726037,0.1527533871307259/  !!! weight coefficient a_k

       c1=(b-a)/2.0
       c2=c1+a
       gl=0.0
       do i=1,n
              c0=c1*t(i)
              t1=c2+c0
              t2=c2-c0
              f1=f(t1)
              f2=f(t2)
              gl=gl+w(i)*(f1+f2)
       end do
       gl=c1*gl

       return
       end subroutine GLI



       function f(x)
       implicit none
       real(8) :: f, x

       f=1.0/(1.0+x) !!! Function for integration

       return
       end function f
