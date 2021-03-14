#define linear_interpolation
!~#define quadratic_interpolation
!~#define cubic_interpolation

#define outputData

    program main 
    implicit none
    real(8), parameter :: Pi=4.0d0*datan(1.0d0)
    integer, parameter :: nx=101
    real(8) :: xExact(nx), uExact(nx)
    
#ifdef linear_interpolation
    integer, parameter :: order=1
#endif
#ifdef quadratic_interpolation
    integer, parameter :: order=2
#endif
#ifdef cubic_interpolation
    integer, parameter :: order=3
#endif

    integer, parameter :: meshX=10
    real(8) :: xMesh(meshX), uMesh(meshX)

    integer, parameter :: particleX=100
    real(8) :: xInterpolated(particleX), uInterpolated(particleX)
    
    integer :: iLoc
    integer :: i, j
    real(8) :: dx, dx2
    integer :: findLoc
    integer :: locate
    real(8) :: errorL1, errorL2
    integer :: errorNum
    character*24 ctime, string
    INTEGER*4  time
    real(kind=8) :: start, finish
    
    !~ write(*,*) "Pi=",Pi
    !~ write(*,*) "nint(Pi+0.5)=", nint(Pi+0.5)
    !~ write(*,*) "int(Pi+0.5)=", int(Pi+0.5)
    
    string = ctime( time() )
    write(*,*) 'Start: ', string
    call CPU_TIME(start)
    
#ifdef linear_interpolation
    write(*,*) "I am linear interpolation!"
#endif
#ifdef quadratic_interpolation
    write(*,*) "I am quadratic interpolation!"
#endif
#ifdef cubic_interpolation
    write(*,*) "I am cubic interpolation!"
#endif
    write(*,*) "    "
    
    !--exact solution
    do i=1,nx
        xExact(i) = dble(i-1)/dble(nx-1)*2.0d0*Pi
        uExact(i) = dsin(xExact(i))
    enddo
#ifdef outputData
    open(unit=01,file="exact.dat",status="unknown")
    do i=1,nx
        write(01,*) xExact(i), uExact(i)
    enddo
    close(01)
#endif
    
    write(*,*) "meshX:", meshX
    !--raw data points
    dx = 2.0d0*Pi/dble(meshX)
    do i=1,meshX
        xMesh(i) = dx/2.0d0+dble(i-1)*dx
        uMesh(i) = dsin(xMesh(i))
    enddo
#ifdef outputData
    open(unit=02,file="raw.dat",status="unknown")
    do i=1,meshX
        write(02,*) xMesh(i), uMesh(i)
    enddo
    close(02)
#endif
    
    write(*,*) "particleX:", particleX
    !--interpolated data points
    dx2 = 2.0d0*Pi/dble(particleX-1)
    errorL1 = 0.0d0
    errorL2 = 0.0d0
    errorNum = 0
    do i=1,particleX
        xInterpolated(i) = dble(i-1)*dx2
        if( (xInterpolated(i).LT.xMesh(1)).OR.(xInterpolated(i).GT.xMesh(meshX)) ) then
            !~ write(*,*) "xInterpolated(i) is out of the range of X1"
            !~ write(*,*) "i =", i
            !~ write(*,*) "xInterpolated(i) =", xInterpolated(i)
            !~ write(*,*) "xMesh_min =", xMesh(1),"  , xMesh_mmeshX = ", xMesh(meshX)
            !~ write(*,*) "    "
            uInterpolated(i) = 0.0d0
        ELSE
#ifdef linear_interpolation               
            iLoc = locate(xMesh, xInterpolated(i), meshX) 
            if( (iLoc.GE.1).AND.(iLoc.LE.meshX-1) ) then
                call LagrangeInterpolation(xMesh(iLoc:iLoc+1), uMesh(iLoc:iLoc+1), xInterpolated(i), uInterpolated(i), order)
            else
                write(*,*) "check boundary, iLoc=", iLoc
                stop
            endif
#endif

#ifdef quadratic_interpolation
            iLoc = locate(xMesh, xInterpolated(i), meshX)             
            if( (iLoc.GE.2).AND.(iLoc.LE.meshX-1) ) then
                call LagrangeInterpolation(xMesh(iLoc-1:iLoc+1), uMesh(iLoc-1:iLoc+1), xInterpolated(i), uInterpolated(i), order)
            elseif( (iLoc.LT.2) ) then
                if(iLoc.LE.0) then
                    write(*,*) "check boundary, iLoc=", iLoc
                    stop
                endif
                call LagrangeInterpolation(xMesh(iLoc:iLoc+2), uMesh(iLoc:iLoc+2), xInterpolated(i), uInterpolated(i), order)
            else
                write(*,*) "check boundary, iLoc=", iLoc
                stop
            endif
#endif

#ifdef cubic_interpolation
            iLoc = locate(xMesh, xInterpolated(i), meshX) 
            if( (iLoc.GE.2).AND.(iLoc.LE.meshX-2) ) then
                call LagrangeInterpolation(xMesh(iLoc-1:iLoc+2), uMesh(iLoc-1:iLoc+2), xInterpolated(i), uInterpolated(i), order)
            elseif( (iLoc.LT.2) ) then
                if(iLoc.LE.0) then
                    write(*,*) "check boundary, iLoc=", iLoc
                    stop
                endif
                call LagrangeInterpolation(xMesh(iLoc:iLoc+3), uMesh(iLoc:iLoc+3), xInterpolated(i), uInterpolated(i), order)
            elseif( (iLoc.GT.meshX-2) ) then
                if(iLoc.GE.meshX) then
                    write(*,*) "check boundary, iLoc=", iLoc
                    stop
                endif
                call LagrangeInterpolation(xMesh(iLoc-2:iLoc+1), uMesh(iLoc-2:iLoc+1), xInterpolated(i), uInterpolated(i), order)
            else
                write(*,*) "check boundary, iLoc=", iLoc
                stop
            endif
#endif

            errorL1= errorL1+abs(uInterpolated(i)-dsin(xInterpolated(i)))
            errorL2 = errorL2+(uInterpolated(i)-dsin(xInterpolated(i)))**2.0d0
            errorNum = errorNum+1
            
        endif
    enddo
    
#ifdef outputData
    open(unit=03,file="interpolated.dat",status="unknown")
    do i=1,particleX
        write(03,*) xInterpolated(i), uInterpolated(i)
    enddo
    close(03)
#endif
    
    call CPU_TIME(finish)
    write(*,*) "Time (CPU) = ", real(finish-start), "s"
    write(*,*) "L1 error=", errorL1/dble(errorNum)
    write(*,*) "L2 error=", dsqrt(errorL2/dble(errorNum))
    
    string = ctime( time() )
    write(*,*) 'End:   ', string
    
    stop
    end program main
    
!--Given an array xx(1:N), and given a value x, returns a value j
!----such that x is between xx(j) and xx(j+1).
!----N = size(xx)
!----xx must be monotonic, either increasing or decreasing.
!----j=0 or j=N is returned to indicate that x is out of range
    function locate(xx, x, N)
    implicit none
    real(8) :: xx(1:N)
    real(8) :: x
    integer :: locate
    integer :: N, jLower, jMedium, jUpper
    logical :: ascnd
    
    ascnd = (xx(n).GE.xx(1))
    jLower = 0
    jUpper = N+1
    do 
        if((jUpper-jLower).LE.1) exit
        jMedium = (jUpper+jLower)/2
        if( ascnd.EQV.(x.GE.xx(jMedium)) ) then
            jLower = jMedium
        else
            jUpper = jMedium
        endif
    enddo
    if (x==xx(1)) then
        locate = 1
    elseif (x == xx(N)) then
        locate = N-1
    else
        locate = jLower
    endif
    
    end function locate
    
    
    subroutine LagrangeInterpolation(pointX, pointU, point0X, point0U, order)
    implicit none
    integer :: order
    real(8) :: pointX(1:order+1), pointU(1:order+1)
    real(8) :: point0X, point0U
    REAL(8) :: TEMP
    integer :: j, k

    point0U = 0.0d0
    do k=0,order
        temp = 1.0d0
        do j=0,order
            if(j.NE.K) temp=temp*(point0X-pointX(j+1))/(pointX(k+1)-pointX(j+1))
        enddo
        point0U= point0U+temp*pointU(k+1)
    enddO
                
    return
    end subroutine LagrangeInterpolation

