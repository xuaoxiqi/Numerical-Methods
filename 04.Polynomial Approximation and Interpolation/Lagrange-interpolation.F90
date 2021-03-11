#define linear_interpolation
!~#define quadratic_interpolation
!~#define cubic_interpolation

#define outputData

    program main 
    implicit none
    real(8), parameter :: Pi=4.0d0*datan(1.0d0)
    integer, parameter :: nx=101
    real(8) :: xExact(nx), yExact(nx)
    
    integer, parameter :: ax=10
    real(8) :: xRaw(ax), yRaw(ax)
    
    integer, parameter :: bx=6
    real(8) :: xInterpolated(bx), yInterpolated(bx)
    
    integer :: i, j
    real(8) :: dx, dx2
    integer :: xLoc, xLoc_plusOne
    real(8) :: yLoc, yLoc_plusOne
    real(8) :: point0X, point0Y
    real(8) :: point1X, point1Y, point2X, point2Y,  point3X, point3Y,  point4X, point4Y
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
        yExact(i) = dsin(xExact(i))
    enddo
#ifdef outputData
    open(unit=01,file="exact.dat",status="unknown")
    do i=1,nx
        write(01,*) xExact(i), yExact(i)
    enddo
    close(01)
#endif
    
    !--raw data points
    dx = 2.0d0*Pi/dble(ax)
    do i=1,ax
        xRaw(i) = dx/2.0d0+dble(i-1)*dx
        yRaw(i) = dsin(xRaw(i))
    enddo
#ifdef outputData
    open(unit=02,file="raw.dat",status="unknown")
    do i=1,ax
        write(02,*) xRaw(i), yRaw(i)
    enddo
    close(02)
#endif
    
    !--interpolated data points
    dx2 = 2.0d0*Pi/dble(bx-1)
    errorL1 = 0.0d0
    errorL2 = 0.0d0
    errorNum = 0
    do i=1,bx
        xInterpolated(i) = dble(i-1)*dx2
        if( (xInterpolated(i).LT.xRaw(1)).OR.(xInterpolated(i).GT.xRaw(ax)) ) then
            !~ write(*,*) "xInterpolated(i) is out of the range of X1"
            !~ write(*,*) "i =", i
            !~ write(*,*) "xInterpolated(i) =", xInterpolated(i)
            !~ write(*,*) "xRaw_min =", xRaw(1),"  , xRaw_max = ", xRaw(ax)
            !~ write(*,*) "    "
            yInterpolated(i) = 0.0d0
        ELSE
#ifdef linear_interpolation               
            j = locate(xRaw, xInterpolated(i), ax) 
            if(j.EQ.ax) then
                call linearInterpolation(xRaw(j-1), yRaw(j-1), xRaw(j), yRaw(j), xInterpolated(i), yInterpolated(i))
            else
                call linearInterpolation(xRaw(j), yRaw(j), xRaw(j+1), yRaw(j+1), xInterpolated(i), yInterpolated(i))
            endif
            errorL1= errorL1+abs(yInterpolated(i)-dsin(xInterpolated(i)))
            errorL2 = errorL2+(yInterpolated(i)-dsin(xInterpolated(i)))**2.0d0
            errorNum = errorNum+1
#endif

#ifdef quadratic_interpolation
            j = locate(xRaw, xInterpolated(i), ax) 
            if(j.EQ.1) then
                call quadraticInterpolation(xRaw(j), yRaw(j), xRaw(j+1), yRaw(j+1), xRaw(j+2), yRaw(j+2), xInterpolated(i), yInterpolated(i))
            elseif(j.EQ.ax) then
                call quadraticInterpolation(xRaw(j-2), yRaw(j-2), xRaw(j-1), yRaw(j-1), xRaw(j), yRaw(j), xInterpolated(i), yInterpolated(i))
            else
                call quadraticInterpolation(xRaw(j-1), yRaw(j-1), xRaw(j), yRaw(j), xRaw(j+1), yRaw(j+1), xInterpolated(i), yInterpolated(i))
            endif
            errorL1= errorL1+abs(yInterpolated(i)-dsin(xInterpolated(i)))
            errorL2 = errorL2+(yInterpolated(i)-dsin(xInterpolated(i)))**2.0d0
            errorNum = errorNum+1
#endif

#ifdef cubic_interpolation
            j = locate(xRaw, xInterpolated(i), ax) 
            if(j.EQ.1) then
                call cubicInterpolation(xRaw(j), yRaw(j), xRaw(j+1), yRaw(j+1),  &
                                            xRaw(j+2), yRaw(j+2), xRaw(j+3), yRaw(j+3), xInterpolated(i), yInterpolated(i))
            elseif(j.EQ.ax) then
                call cubicInterpolation(xRaw(j-3), yRaw(j-3), xRaw(j-2), yRaw(j-2),  &
                                            xRaw(j-1), yRaw(j-1), xRaw(j), yRaw(j), xInterpolated(i), yInterpolated(i))
            elseif(j.EQ.ax-1) then
                call cubicInterpolation(xRaw(j-2), yRaw(j-2), xRaw(j-1), yRaw(j-1),  &
                                            xRaw(j), yRaw(j), xRaw(j+1), yRaw(j+1), xInterpolated(i), yInterpolated(i))
            else
                call cubicInterpolation(xRaw(j-1), yRaw(j-1), xRaw(j), yRaw(j),  &
                                            xRaw(j+1), yRaw(j+1), xRaw(j+2), yRaw(j+2), xInterpolated(i), yInterpolated(i))
            endif
            errorL1= errorL1+abs(yInterpolated(i)-dsin(xInterpolated(i)))
            errorL2 = errorL2+(yInterpolated(i)-dsin(xInterpolated(i)))**2.0d0
            errorNum = errorNum+1
#endif
        endif
    enddo
    
#ifdef outputData
    open(unit=03,file="interpolated.dat",status="unknown")
    do i=1,bx
        write(03,*) xInterpolated(i), yInterpolated(i)
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
        if( ascnd.EQ.(x.GE.xx(jMedium)) ) then
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
    
    
    subroutine linearInterpolation(point1X, point1Y, point2X, point2Y, point0X, point0Y)
    implicit none
    real(8) :: point1X, point1Y, point2X, point2Y
    real(8) :: point0X, point0Y
    
    !--Straightforward implementation
    point0Y = (point0X-point2X)/(point1X-point2X)*point1Y+(point0X-point1X)/(point2X-point1X)*point2Y
    
    return
    end subroutine linearInterpolation
    
    
    subroutine quadraticInterpolation(point1X, point1Y, point2X, point2Y, point3X, point3Y, point0X, point0Y)
    implicit none
    real(8) :: point1X, point1Y, point2X, point2Y, point3X, point3Y
    real(8) :: point0X, point0Y
    
    point0Y = (point0X-point2X)*(point0X-point3X)/(point1X-point2X)/(point1X-point3X)*point1Y &
                +(point0X-point1X)*(point0X-point3X)/(point2X-point1X)/(point2X-point3X)*point2Y &
                +(point0X-point1X)*(point0X-point2X)/(point3X-point1X)/(point3X-point2X)*point3Y 
                
    return
    end subroutine quadraticInterpolation
    

    subroutine cubicInterpolation(point1X, point1Y, point2X, point2Y, point3X, point3Y, point4X, point4Y, point0X, point0Y)
    implicit none
    real(8) :: point1X, point1Y, point2X, point2Y, point3X, point3Y, point4X, point4Y
    real(8) :: point0X, point0Y
    
    point0Y = (point0X-point2X)*(point0X-point3X)*(point0X-point4X) &
                    /(point1X-point2X)/(point1X-point3X)/(point1X-point4X)*point1Y &
                +(point0X-point1X)*(point0X-point3X)*(point0X-point4X) &
                    /(point2X-point1X)/(point2X-point3X)/(point2X-point4X)*point2Y &
                +(point0X-point1X)*(point0X-point2X)*(point0X-point4X) &
                    /(point3X-point1X)/(point3X-point2X)/(point3X-point4X)*point3Y &
                +(point0X-point1X)*(point0X-point2X)*(point0X-point3X) &
                    /(point4X-point1X)/(point4X-point2X)/(point4X-point3X)*point4Y 
                
    return
    end subroutine cubicInterpolation

