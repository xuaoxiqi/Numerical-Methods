!~#define linear_interpolation
!~#define quadratic_interpolation
#define cubic_interpolation

#define outputData

    program main 
    implicit none
    real(8), parameter :: Pi=4.0d0*datan(1.0d0)
    integer, parameter :: nx=101, ny=81
    real(8) :: xExact(nx), yExact(ny), uExact(nx,ny)
    
    integer, parameter :: meshX=50, meshY=40
    real(8) :: xMesh(meshX), yMesh(meshY), uMesh(meshX,meshY)
    
    integer, parameter :: particleX=200, particleY=160
    real(8) :: xInterpolated(particleX), yInterpolated(particleY), uInterpolated(particleX,particleY)
    
    integer :: i, j
    integer :: iLoc, jLoc
    real(8) :: dx, dy, dx2, dy2
    real(8) :: point0X, point0U
    real(8) :: point1X, point1U, point2X, point2U,  point3X, point3U,  point4X, point4U
    integer :: findLoc
    integer :: locate
    
    real(8) :: errorL1, errorL2
    integer :: errorNum
    
    character*24 ctime, string
    INTEGER*4  time
    real(kind=8) :: start, finish
    real(8) :: tempY(1:4)
    
    
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
    enddo
    do j=1,ny
        yExact(j) = dble(j-1)/dble(ny-1)*2.0d0*Pi
    enddo
    do j=1,ny
        do i=1,nx
            uExact(i,j) = dsin(xExact(i))*dcos(yExact(j))
        enddo
    enddo
#ifdef outputData
    open(unit=01,file="exact.dat",status="unknown")
    write(01,*) 'TITLE="Lid Driven Cavity"'
    write(01,*) 'VARIABLES="X" "Y" "U" '
    write(01,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(01,100) xExact(i), yExact(j), uExact(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format(' ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(01)
#endif
    
    !--raw data points
    dx = 2.0d0*Pi/dble(meshX)
    dy = 2.0d0*Pi/dble(meshY)
    do i=1,meshX
        xMesh(i) = dx/2.0d0+dble(i-1)*dx
    enddo
    do j=1,meshY
        yMesh(j) = dy/2.0d0+dble(j-1)*dy
    enddo
    do j=1,meshY
        do i=1,meshX
            uMesh(i,j) = dsin(xMesh(i))*dcos(yMesh(j))
        enddo
    enddo
#ifdef outputData
    open(unit=02,file="raw.dat",status="unknown")
    write(02,*) 'TITLE="Lid Driven Cavity"'
    write(02,*) 'VARIABLES="X" "Y" "uMesh" '
    write(02,201) meshX, meshY
    do j=1,meshY
        do i=1,meshX
            write(02,200) xMesh(i), yMesh(j), uMesh(i,j)
        enddo
    enddo
200 format(1x,2(e11.4,' '),10(e13.6,' '))
201 format(' ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)
#endif
    
    !--interpolated data points
    dx2 = 2.0d0*Pi/dble(particleX-1)
    dy2 = 2.0d0*Pi/dble(particleY-1)
    errorL1 = 0.0d0
    errorL2 = 0.0d0
    errorNum = 0
    do i=1,particleX
        xInterpolated(i) = dble(i-1)*dx2
    enddo
    do j=1,particleY
        yInterpolated(j)= dble(j-1)*dy2
    enddo
    do j=1,particleY
        do i=1,particleX
            if( (xInterpolated(i).LT.xMesh(1)).OR.(xInterpolated(i).GT.xMesh(meshX)).OR.(yInterpolated(j).LT.yMesh(1)).OR.(yInterpolated(j).GT.yMesh(meshY)) ) then
                !~ write(*,*) "xInterpolated(i) is out of the range of X1"
                !~ write(*,*) "i =", i
                !~ write(*,*) "xInterpolated(i) =", xInterpolated(i)
                !~ write(*,*) "xMesh_min =", xMesh(1),"  , xMesh_mmeshX = ", xMesh(meshX)
                !~ write(*,*) "    "
                !~uInterpolated(i,j) = 0.0d0
                !~write(*,*) "i,j=", i, j
            ELSE
#ifdef linear_interpolation
                iLoc = locate(xMesh, xInterpolated(i), meshX)
                jLoc = locate(yMesh, yInterpolated(j), meshY) 
                if( (iLoc.EQ.meshX).OR.(jLoc.EQ.meshY) ) then
                    write(*,*) "check (iLoc.EQ.meshX).OR.(jLoc.EQ.meshY)"
                    stop
                else
                    call linearInterpolation(yMesh(jLoc), uMesh(iLoc,jLoc), yMesh(jLoc+1), uMesh(iLoc,jLoc+1), yInterpolated(j), tempY1)
                    call linearInterpolation(yMesh(jLoc), uMesh(iLoc+1,jLoc), yMesh(jLoc+1), uMesh(iLoc+1,jLoc+1), yInterpolated(j), tempY2)
                    call linearInterpolation(xMesh(iLoc), tempY1, xMesh(iLoc+1), tempY2, xInterpolated(i), uInterpolated(i,j))
                endif
                
                errorL1= errorL1+abs(uInterpolated(i,j)-dsin(xInterpolated(i))*dcos(yInterpolated(j)))
                errorL2 = errorL2+(uInterpolated(i,j)-dsin(xInterpolated(i))*dcos(yInterpolated(j)))**2.0d0
                errorNum = errorNum+1
#endif

#ifdef cubic_interpolation
                iLoc = locate(xMesh, xInterpolated(i), meshX)
                jLoc = locate(yMesh, yInterpolated(j), meshY) 
                if( (iLoc.GE.2).AND.(jLoc.GE.2).AND.(iLoc.LE.meshX-2).AND.(jLoc.LE.meshY-2) ) then
                    !--Bulk
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc-1, jLoc-1:jLoc+2), yInterpolated(j), tempY(1))
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc,    jLoc-1:jLoc+2), yInterpolated(j), tempY(2))
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc+1,jLoc-1:jLoc+2), yInterpolated(j), tempY(3))
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc+2,jLoc-1:jLoc+2), yInterpolated(j), tempY(4))
                    
                    call cubicInterpolation(xMesh(iLoc-1:iLoc+2), tempY(1:4), xInterpolated(i), uInterpolated(i,j))

                elseif( (iLoc.LT.2).AND.(jLoc.GE.2).AND.(jLoc.LE.meshY-2) ) then
                    !----Left wall
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc,    jLoc-1:jLoc+2), yInterpolated(j), tempY(1))
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc+1,jLoc-1:jLoc+2), yInterpolated(j), tempY(2))
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc+2,jLoc-1:jLoc+2), yInterpolated(j), tempY(3))
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc+3,jLoc-1:jLoc+2), yInterpolated(j), tempY(4))
                    
                    call cubicInterpolation(xMesh(iLoc:iLoc+3), tempY(1:4), xInterpolated(i), uInterpolated(i,j))
                
                elseif( (iLoc.GE.2).AND.(jLoc.LT.2).AND.(iLoc.LE.meshX-2) ) then
                    !----Bottom wall
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc-1, jLoc:jLoc+3), yInterpolated(j), tempY(1))
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc,    jLoc:jLoc+3), yInterpolated(j), tempY(2))
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc+1,jLoc:jLoc+3), yInterpolated(j), tempY(3))
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc+2,jLoc:jLoc+3), yInterpolated(j), tempY(4))
                    
                    call cubicInterpolation(xMesh(iLoc-1:iLoc+2), tempY(1:4), xInterpolated(i), uInterpolated(i,j))

                elseif( (jLoc.GE.2).AND.(iLoc.GT.meshX-2).AND.(jLoc.LE.meshY-2) ) then
                    !----Right wall
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc-2, jLoc-1:jLoc+2), yInterpolated(j), tempY(1))
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc-1, jLoc-1:jLoc+2), yInterpolated(j), tempY(2))
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc,    jLoc-1:jLoc+2), yInterpolated(j), tempY(3))
                    call cubicInterpolation(yMesh(jLoc-1:jLoc+2), uMesh(iLoc+1,jLoc-1:jLoc+2), yInterpolated(j), tempY(4))
                    
                    call cubicInterpolation(xMesh(iLoc-2:iLoc+1), tempY(1:4), xInterpolated(i), uInterpolated(i,j))
                
                elseif( (iLoc.GE.2).AND.(iLoc.LE.meshX-2).AND.(jLoc.GT.meshY-2) ) then
                    !----Top wall
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc-1, jLoc-2:jLoc+1), yInterpolated(j), tempY(1))
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc,    jLoc-2:jLoc+1), yInterpolated(j), tempY(2))
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc+1,jLoc-2:jLoc+1), yInterpolated(j), tempY(3))
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc+2,jLoc-2:jLoc+1), yInterpolated(j), tempY(4))
                
                    call cubicInterpolation(xMesh(iLoc-1:iLoc+2), tempY(1:4), xInterpolated(i), uInterpolated(i,j))
                
                elseif( (iLoc.LT.2).AND.(jLoc.LT.2) ) then
                    !------Left&Bottom corner
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc,    jLoc:jLoc+3), yInterpolated(j), tempY(1))
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc+1,jLoc:jLoc+3), yInterpolated(j), tempY(2))
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc+2,jLoc:jLoc+3), yInterpolated(j), tempY(3))
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc+3,jLoc:jLoc+3), yInterpolated(j), tempY(4))
                    
                    call cubicInterpolation(xMesh(iLoc:iLoc+3), tempY(1:4), xInterpolated(i), uInterpolated(i,j))
                
                elseif( (iLoc.LT.2).AND.(jLoc.GT.meshY-2) ) then
                    !------Left&Top corner
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc,    jLoc-2:jLoc+1), yInterpolated(j), tempY(1))
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc+1,jLoc-2:jLoc+1), yInterpolated(j), tempY(2))
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc+2,jLoc-2:jLoc+1), yInterpolated(j), tempY(3))
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc+3,jLoc-2:jLoc+1), yInterpolated(j), tempY(4))
                    
                    call cubicInterpolation(xMesh(iLoc:iLoc+3), tempY(1:4), xInterpolated(i), uInterpolated(i,j))
                
                elseif( (jLoc.LT.2).AND.(iLoc.GT.meshX-2) ) then
                    !------Right&Bottom corner
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc-2, jLoc:jLoc+3), yInterpolated(j), tempY(1))
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc-1, jLoc:jLoc+3), yInterpolated(j), tempY(2))
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc,    jLoc:jLoc+3), yInterpolated(j), tempY(3))
                    call cubicInterpolation(yMesh(jLoc:jLoc+3), uMesh(iLoc+1,jLoc:jLoc+3), yInterpolated(j), tempY(4))
                    
                    call cubicInterpolation(xMesh(iLoc-2:iLoc+1), tempY(1:4), xInterpolated(i), uInterpolated(i,j))
                
                elseif( (iLoc.GT.meshX-2).AND.(jLoc.GT.meshY-2) ) then
                    !------Right&Top corner
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc-2, jLoc-2:jLoc+1), yInterpolated(j), tempY(1))
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc-1, jLoc-2:jLoc+1), yInterpolated(j), tempY(2))
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc,    jLoc-2:jLoc+1), yInterpolated(j), tempY(3))
                    call cubicInterpolation(yMesh(jLoc-2:jLoc+1), uMesh(iLoc+1,jLoc-2:jLoc+1), yInterpolated(j), tempY(4))
                    
                    call cubicInterpolation(xMesh(iLoc-2:iLoc+1), tempY(1:4), xInterpolated(i), uInterpolated(i,j))
                
                else
                    write(*,*) "Check boundary..."
                    write(*,*) "iLoc=,", iLoc, ", jLoc=", jLoc
                    stop
                endif
                
                errorL1= errorL1+abs(uInterpolated(i,j)-dsin(xInterpolated(i))*dcos(yInterpolated(j)))
                errorL2 = errorL2+(uInterpolated(i,j)-dsin(xInterpolated(i))*dcos(yInterpolated(j)))**2.0d0
                errorNum = errorNum+1
#endif
            endif
        enddo
    enddo
    
#ifdef outputData  
    open(unit=03,file="interpolated.dat",status="unknown")
    write(03,*) 'TITLE="Lid Driven Cavity"'
    write(03,*) 'VARIABLES="X" "Y" "uInterpolated" '
    write(03,301) particleX, particleY
    do j=1,particleY
        do i=1,particleX
            write(03,300) xInterpolated(i), yInterpolated(j), uInterpolated(i,j)
        enddo
    enddo
300 format(1x,2(e11.4,' '),10(e13.6,' '))
301 format(' ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
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
    
!--Given an arrmeshY xx(1:N), and given a value x, returns a value j
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

