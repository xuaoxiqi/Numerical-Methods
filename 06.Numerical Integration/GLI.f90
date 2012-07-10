
!!!    This program integrates function f(x) with respect to x from a to b
!!!    by using Gauss-Legendre numerical integration method(n=20)
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

       program main
       implicit none
       real(8) :: f, gl, a, b
       external f

       a = 0.0d0
       b = 1.0d0
       write(*,*) 'Gauss-Legendre Method(N=20):'
       write(*,*) 'f(x) = 1.0/(1.0+x)'
       write(*,*) 'Integration on [',a,b,']'
       call GLI(a,b,gl)
       write(*,*) 'Result:',gl

       stop
       end program main



       subroutine GLI(a,b,gl)
       implicit none
       integer, parameter :: k=10
       integer :: i
       real(8) :: x_k(k),A_k(k)
       real(8) :: a, b, c0, c1, c2, t1, t2, f, f1, f2, gl

!!!     n = 20, k =10
       data x_k/0.9931285991850949d0,0.9639719272779138d0,0.9122344282513259d0,0.8391169718222188d0,0.7463319064601508d0, &
                0.6360536807265150d0,0.5108670019508271d0,0.3737060887154196d0,0.2277858511416451d0,0.07652652113349734d0/
       data A_k/0.01761400713915212d0,0.04060142980038694d0,0.06267204833410906d0,0.08327674157670475d0,0.1019301198172404d0, &
                0.1181945319615184d0,0.1316886384491766d0,0.1420961093183821d0,0.1491729864726037d0,0.1527533871307259d0/

!!!    n = 10, k=5
!!!    data x_k/0.9739065285171717d0,0.8650633666889845d0,0.6794095682990244d0,0.4333953941292472d0,0.1488743389816312d0/
!!!    data A_k/0.06667134430868814d0,0.1494513491505806d0,0.2190863625159820d0,0.2692667193099964d0,0.2955242247147529d0/

!!!    n = 5, k = 3
!!!    data x_k/0.000000000000000d0,0.538469310105683d0,0.906179845938664d0/
!!!    data A_k/0.568888888888889d0,0.478628670499366d0,0.236926885056189d0/

!!!    n = 4, k = 2
!!!    data x_k/0.339981043584856d0,0.861136311594053d0/
!!!    data A_k/0.347854845137454d0,0.652145154862546d0/

!!!    n = 3, k = 2
!!!    data x_k/0.000000000000000d0,0.774596669241483d0/
!!!    data A_k/0.888888888888889d0,0.555555555555556d0/

!!!    n = 2, k = 1
!!!    data x_k/0.577350269189626d0/
!!!    data A_k/1.000000000000000d0/

       c1=(b-a)/2.0d0
       c2=c1+a
       gl=0.0d0
       do i=1,k
              c0=c1*x_k(i)
              t1=c2+c0
              t2=c2-c0
              f1=f(t1)
              f2=f(t2)
              gl=gl+A_k(i)*(f1+f2)
       end do
       gl=c1*gl

       return
       end subroutine GLI



       function f(x)
       implicit none
       real(8) :: f, x

       f=1.0d0/(1.0d0+x) !!! Function for integration

       return
       end function f
