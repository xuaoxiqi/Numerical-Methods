
!!!    This program finds eigenvalues and eigenvectors of a real symmetric matrix by using Jacobi Method.
!!!    This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License. 
!!!    Ao Xu, Profiles: <http://www.linkedin.com/pub/ao-xu/30/a72/a29>

    program main
    implicit none
    integer, parameter :: n=4
    integer i, j
    integer :: nRot ! the numer of Jacobi rotations that were required
    real, parameter:: eps=1.0e-09
    real(8) :: a(n,n) ! a is a real symmetric matrix 
    real(8) :: v(n,n) ! v's colums contain the normalized eigenvectors of a
    real(8) :: d(n) ! eigenvalues of a in its first n elements
    
    a(1,:)=(/1.0, 2.0, 1.0, 2.0/)
    a(2,:)=(/2.0, 2.0, -1.0, 1.0/)
    a(3,:)=(/1.0, -1.0, 1.0, 1.0/)
    a(4,:)=(/2.0, 1.0, 1.0, 1.0/)

    write(*,*) 'Jacobi Method:'
    write(*,*) 'All the elements of matrix A:'
    write(*,'(4f9.4)') ((a(i,j),j=1,n),i=1,n)
    call Jacobi(a,4,d,v,nRot)
    write(*,*) "nRot=", nRot
    
    do i=1,n
        write(*,*) '*************************************'
        write(*,*)
        write (*,*) 'Eigenvalues:'
        write (*,'(4f9.4)') d(i)
        write(*,*) 'Eigenvector:'
        write(*,15) (v(j,i),j=1,n)
15     format(1x,'[',f9.6,',',f9.6,',',f9.6,',',f9.6,']')
    enddo

    stop
    end program main



    subroutine Jacobi(a,n,d,v,nRot)
    implicit none
    integer :: n, nRot
    real(8) :: a(n,n), d(n), v(n,n)
    integer :: i, ip, iq, j
    real(8) :: c, g, h, s, sm, t, tau, theta, tresh
    real(8) :: b(n), z(n)
    
    ! initialize v  to the identity matrix
    do ip=1,n
        do iq=1,n
            v(ip,iq) = 0.0d0
        enddo
        v(ip,ip) = 1.0d0
    enddo
    
    ! initizlize b and d to the diagonal of a
    do ip=1,n
        b(ip) = a(ip,ip)
        d(ip) = b(ip)
        z(ip) = 0.0d0
    enddo
    
    nRot = 0
    
    do i=1,50
        ! sum off-diagonal elements
        sm = 0.0d0
        do ip=1,n-1
            do iq=ip+1,n
                sm = sm+DABS(a(ip,iq))
            enddo
        enddo
        ! Normal return ...
        if(sm.EQ.0.0d0) return
        ! On the first three sweeps, rotate only if tresh exceeded.
        if(i.LT.4) then
            tresh = 0.2d0*sm/n**2
        else
            tresh = 0.0d0
        endif
        
        do ip=1,n-1
            do iq=ip+1,n
                g = 100.0d0*DABS(a(ip,iq))
                ! After four sweeps, skip the rotation if the off-diagnoal element is small
                if( (i.GT.4).AND.(ABS(d(ip))+g.EQ.ABS(d(ip))).AND.(ABS(d(iq))+g.eq.ABS(d(iq))) ) then
                    a(ip,iq) = 0.0d0
                elseif(ABS(a(ip,iq)).GT.tresh) then
                    h = d(iq)-d(ip)
                    if(ABS(h)+g.EQ.ABS(h)) then
                        t = a(ip,iq)/h  ! t=1/(2theta)
                    else
                        theta = 0.5d0*h/a(ip,iq)
                        t = 1.0d0/(ABS(theta)+SQRT(1.0d0+theta**2))
                        if(theta.LT.0.0d0) t = -t
                    endif
                    c = 1.0d0/SQRT(1.0d0+t**2)
                    s = t*c
                    tau = s/(1.0d0+c)
                    h = t*a(ip,iq)
                    z(ip) = z(ip)-h
                    z(iq) = z(iq)+h
                    d(ip) = d(ip)-h
                    d(iq) = d(iq)+h
                    a(ip,iq) = 0.0d0
                    
                    ! case of ratations 1<j<p
                    do j=1,ip-1
                        g = a(j,ip)
                        h = a(j,iq)
                        a(j,ip) = g-s*(h+g*tau)
                        a(j,iq) = h+s*(g-h*tau)
                    enddo
                    
                    ! case of ratations p<j<q
                    do j=ip+1,iq-1
                        g = a(ip,j)
                        h = a(j,iq)
                        a(ip,j) = g-s*(h+g*tau)
                        a(j,iq) = h+s*(g-h*tau)
                    enddo
                    
                    ! case of rotations q<j<n
                    do j=iq+1,n
                        g = a(ip,j)
                        h = a(iq,j)
                        a(ip,j) = g-s*(h+g*tau)
                        a(iq,j) = h+s*(g-h*tau)
                    enddo
                    
                    do j=1,n
                        g = v(j,ip)
                        h = v(j,iq)
                        v(j,ip) = g-s*(h+g*tau)
                        v(j,iq) = h+s*(g-h*tau)
                    enddo
                    
                    nRot = nRot+1
                endif
            enddo
        enddo
        
        do ip=1,n
            b(ip) = b(ip)+z(ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0 ! reinitizlize z
        enddo
    enddo
                        
    return
    end subroutine Jacobi
