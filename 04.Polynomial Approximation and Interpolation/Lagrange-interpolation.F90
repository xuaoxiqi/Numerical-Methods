!~#define linear_interpolation
!~#define quadratic_interpolation
#define cubic_interpolation

#define outputData

    program main 
    implicit none
    real(8), parameter :: Pi=4.0d0*datan(1.0d0)
    integer, parameter :: nx=101
    real(8) :: xExact(nx), uExact(nx)
    
    integer, parameter :: meshX=10
    real(8) :: xMesh(meshX), yMesh(meshX)
    
    integer, parameter :: particleX=100
    real(8) :: xInterpolated(particleX), uInterpolated(particleX)
    integer :: iLoc
    integer :: i, j
    real(8) :: dx, dx2
    real(8) :: point0X, point0U
    real(8) :: point1X, point1U, point2X, point2U,  point3X, point3U,  point4X, point4U
    integer :: findLoc
    integer :: locate
    
    real(8) :: errorL1, errorL2
    integer :: errorNum
    
    character*24 ctime, string
    INTEGER*4  time
    real(kind=8) :: start, finish
    
    
    write(*,*) "Pi=",Pi
    write(*,*) "nint(Pi+0.5)=", nint(Pi+0.5)
    write(*,*) "int(Pi+0.5)=", int(Pi+0.5)
    
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
    
    !--raw data points
    dx = 2.0d0*Pi/dble(meshX)
    do i=1,meshX
        xMesh(i) = dx/2.0d0+dble(i-1)*dx
        yMesh(i) = dsin(xMesh(i))
    enddo
#ifdef outputData
    open(unit=02,file="raw.dat",status="unknown")
    do i=1,meshX
        write(02,*) xMesh(i), yMesh(i)
    enddo
    close(02)
#endif
    
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
            j = locate(xMesh, xInterpolated(i), meshX) 
            if(j.EQ.meshX) then
                call linearInterpolation(xMesh(j-1), yMesh(j-1), xMesh(j), yMesh(j), xInterpolated(i), uInterpolated(i))
            else
                call linearInterpolation(xMesh(j), yMesh(j), xMesh(j+1), yMesh(j+1), xInterpolated(i), uInterpolated(i))
            endif
            errorL1= errorL1+abs(uInterpolated(i)-dsin(xInterpolated(i)))
            errorL2 = errorL2+(uInterpolated(i)-dsin(xInterpolated(i)))**2.0d0
            errorNum = errorNum+1
#endif

#ifdef quadratic_interpolation
            j = locate(xMesh, xInterpolated(i), meshX) 
            if(j.EQ.1) then
                call quadraticInterpolation(xMesh(j), yMesh(j), xMesh(j+1), yMesh(j+1), xMesh(j+2), yMesh(j+2), xInterpolated(i), uInterpolated(i))
            elseif(j.EQ.meshX) then
                call quadraticInterpolation(xMesh(j-2), yMesh(j-2), xMesh(j-1), yMesh(j-1), xMesh(j), yMesh(j), xInterpolated(i), uInterpolated(i))
            else
                call quadraticInterpolation(xMesh(j-1), yMesh(j-1), xMesh(j), yMesh(j), xMesh(j+1), yMesh(j+1), xInterpolated(i), uInterpolated(i))
            endif
            errorL1= errorL1+abs(uInterpolated(i)-dsin(xInterpolated(i)))
            errorL2 = errorL2+(uInterpolated(i)-dsin(xInterpolated(i)))**2.0d0
            errorNum = errorNum+1
#endif

#ifdef cubic_interpolation
            iLoc = locate(xMesh, xInterpolated(i), meshX) 
            if((iLoc.GE.2).AND.(iLoc.LE.meshX-2)) then
                call cubicInterpolation(xMesh(iLoc-1:iLoc+2), yMesh(iLoc-1:iLoc+2), xInterpolated(i), uInterpolated(i))
            elseif( (iLoc.LT.2) ) then
                if(iLoc.LE.0) then
                    write(*,*) "check boundary, iLoc=", iLoc
                    stop
                endif
                call cubicInterpolation(xMesh(iLoc:iLoc+3), yMesh(iLoc:iLoc+3), xInterpolated(i), uInterpolated(i))
            elseif( (iLoc.GT.meshX-2) ) then
                if(iLoc.GE.meshX) then
                    write(*,*) "check boundary, iLoc=", iLoc
                    stop
                endif
                call cubicInterpolation(xMesh(iLoc-2:iLoc+1), yMesh(iLoc-2:iLoc+1), xInterpolated(i), uInterpolated(i))
            else
                write(*,*) "check boundary, iLoc=", iLoc
                stop
            endif
            errorL1= errorL1+abs(uInterpolated(i)-dsin(xInterpolated(i)))
            errorL2 = errorL2+(uInterpolated(i)-dsin(xInterpolated(i)))**2.0d0
            errorNum = errorNum+1
#endif
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
    
    
    subroutine linearInterpolation(point1X, point1U, point2X, point2U, point0X, point0U)
    implicit none
    real(8) :: point1X, point1U, point2X, point2U
    real(8) :: point0X, point0U
    
    !--Straightforward implementation
    point0U = (point0X-point2X)/(point1X-point2X)*point1U+(point0X-point1X)/(point2X-point1X)*point2U
    
    return
    end subroutine linearInterpolation
    
    
    subroutine quadraticInterpolation(point1X, point1U, point2X, point2U, point3X, point3U, point0X, point0U)
    implicit none
    real(8) :: point1X, point1U, point2X, point2U, point3X, point3U
    real(8) :: point0X, point0U
    
    point0U = (point0X-point2X)*(point0X-point3X)/(point1X-point2X)/(point1X-point3X)*point1U &
                +(point0X-point1X)*(point0X-point3X)/(point2X-point1X)/(point2X-point3X)*point2U &
                +(point0X-point1X)*(point0X-point2X)/(point3X-point1X)/(point3X-point2X)*point3U 
                
    return
    end subroutine quadraticInterpolation
    

    subroutine cubicInterpolation(pointX, pointY, point0X, point0U)
    implicit none
    real(8) :: pointX(1:4)
    real(8) :: pointY(1:4)
    real(8) :: point0X, point0U
    real(8) :: point1X, point1U, point2X, point2U, point3X, point3U, point4X, point4U
    
    point1X = pointX(1)
    point2X = pointX(2)
    point3X = pointX(3)
    point4X = pointX(4)
    
    point1U = pointY(1)
    point2U = pointY(2)
    point3U = pointY(3)
    point4U = pointY(4)
    
    point0U = (point0X-point2X)*(point0X-point3X)*(point0X-point4X) &
                    /(point1X-point2X)/(point1X-point3X)/(point1X-point4X)*point1U &
                +(point0X-point1X)*(point0X-point3X)*(point0X-point4X) &
                    /(point2X-point1X)/(point2X-point3X)/(point2X-point4X)*point2U &
                +(point0X-point1X)*(point0X-point2X)*(point0X-point4X) &
                    /(point3X-point1X)/(point3X-point2X)/(point3X-point4X)*point3U &
                +(point0X-point1X)*(point0X-point2X)*(point0X-point3X) &
                    /(point4X-point1X)/(point4X-point2X)/(point4X-point3X)*point4U 
                
    return
    end subroutine cubicInterpolation

