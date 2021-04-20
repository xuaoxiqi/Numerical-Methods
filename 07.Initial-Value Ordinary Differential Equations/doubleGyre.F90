    module commondata
    implicit none
    real(8), parameter :: Pi=4.0d0*datan(1.0d0)
    integer(8), parameter :: nx=201, ny=101
    real(8), parameter :: tMax=100.0d0
    real(8), parameter :: dt = 0.001d0
    integer(8), parameter :: maxTimeStep=tMax/dt
    real(8), allocatable :: xGrid(:), yGrid(:)
    real(8), allocatable :: u(:,:), v(:,:), timeLocal(:)
    real(8), allocatable :: uParticle(:), vParticle(:)
    real(8), allocatable :: xParticle(:), yParticle(:)
    real(8), parameter :: omega=2.0d0*Pi/10.0d0, myEpsilon=0.1d0, A=0.1d0
    integer(8), parameter :: calRelError=1
    end module commondata
    
!~#define firstEuler
!~#define secondEuler
#define RK4

    program main
    use commondata
    implicit none
    integer(8) :: i, j, k
    real(8) :: f, dfdx
    character(len=100) :: filename
    real(8) :: errorL1, errorL2, errorTime
    integer(8) :: scaleFactor
    integer(8),parameter :: tempMax=tMax/0.00001d0
    real(8), allocatable :: testX(:), testY(:), testTime(:)
    
    allocate (xGrid(nx))
    allocate (yGrid(ny))
    allocate (u(nx,ny))
    allocate (v(nx,ny))
    allocate (timeLocal(0:maxTimeStep))
    allocate (uParticle(0:maxTimeStep))
    allocate (vParticle(0:maxTimeStep))
    allocate (xParticle(0:maxTimeStep))
    allocate (yParticle(0:maxTimeStep))
    allocate (testX(0:tempMax))
    allocate (testY(0:tempMax))
    allocate (testTime(0:tempMax))
    
    write(*,*) "Pi=", Pi
    
    xGrid = 0.0d0
    yGrid = 0.0d0
    do i=1,nx
        xGrid(i) = 2.0d0/(nx-1)*(i-1)
    enddo
    do j=1,ny
        yGrid(j) = 1.0d0/(ny-1)*(j-1)
    enddo

    write(*,*) "dt=", dt
    write(*,*) "tMax=", tMax
    write(filename,*) real(dt)
    filename = adjustl(filename)
    write(*,*) "maxTimeStep=tMax/dt", maxTimeStep
    write(*,*) "        "
    
    !~!~pause
        
#ifdef firstEuler
    write(*,*) "I am firstEuler"
#endif

#ifdef secondEuler
    write(*,*) "I am secondEuler"
#endif

#ifdef RK4
    write(*,*) "I am RK4"
#endif

    timeLocal = 0.0d0
    do k=1,maxTimeStep
        timeLocal(k) = dble(k)*dt
    enddo
    
    xParticle = 0.0d0
    yParticle = 0.0d0
    uParticle = 0.0d0
    vParticle = 0.0d0
    xParticle(0) = 0.5d0
    yParticle(0) = 0.25d0
    do k=0, maxTimeStep-1
        if(MOD(k,2000).EQ.0) write(*,*) "k=",k, "; maxTimeStep=",maxTimeStep, real(real(k)/real(maxTimeStep)*100.0),"%"

            !~ do j=1,ny
                !~ do i=1,nx
                    !~ u(i,j) = -Pi*A*dsin(Pi*f(xGrid(i), time(k)))*dcos(Pi*yGrid(j))
                    !~ v(i,j) = Pi*A*dcos(Pi*f(xGrid(i), time(k)))*dsin(Pi*yGrid(j))*dfdx(xGrid(i), time(k))
                !~ enddo
            !~ enddo
            
#ifdef firstEuler
        call firstEulerUpdate(k)
#endif

#ifdef secondEuler
        call secondEulerUpdate(k)
#endif

#ifdef RK4
        call RK4Update(k)
#endif
    ENDDO
    
    DO K=0,maxTimeStep
        open(unit=02,file='particle'//trim(filename)//'.dat', status="unknown", position="append")
        write(02,*) xParticle(k), yParticle(k), timeLocal(k)
        close(02)
    enddo
    
    !~!~write(*,*) "timeLocal(1)", timeLocal(1)
    
    if(calRelError.EQ.1) then
        scaleFactor = tempMax/maxTimeStep
        write(*,*) "    "
        write(*,*) "tempMax", tempMax
        write(*,*) "scaleFactor=tempMax/timeStepMax", scaleFactor
        !-------------calculate error----------------------
        testX = 0.0d0
        testY = 0.0d0
        testTime = 0.0d0
        open(unit=01, file='benchmark.dat',status="old")
        do i=0,tempMax
            read(01,*) testX(i), testY(i), testTime(i)
        enddo
        close(01)
        
        errorTime = 0.0d0
        errorL1 = 0.0d0
        errorL2 = 0.0d0
        do k=0,maxTimeStep
            errorTime = errorTime+dabs(timeLocal(k)-testTime(k*scaleFactor))
            errorL1 = errorL1+dabs(xParticle(k)-testX(k*scaleFactor))+dabs(yParticle(k)-testY(k*scaleFactor))
            errorL2 = errorL2+(xParticle(k)-testX(k*scaleFactor))**2.0d0 &
                                     +(yParticle(k)-testY(k*scaleFactor))**2.0d0
            !~ write(*,*) "k",k,"; k*scaleFactor",k*scaleFactor,"; errorTime", errorTime
            !~ write(*,*) "timeLocal(k)", timeLocal(k), ";   testTime(k*scaleFactor)",testTime(k*scaleFactor)
            !~ write(*,*) "    "
        enddo
        errorTime = errorTime/maxTimeStep
        errorL1 = errorL1/maxTimeStep
        errorL2 = dsqrt(errorL2/maxTimeStep)
        write(*,*) "errorTime", errorTime
        write(*,*) "errorL1=", errorL1
        write(*,*) "errorL2=", errorL2
    endif
       
    !~!~call output_ASCII()
    
    stop
    end program
    
    
    subroutine firstEulerUpdate(k)
    use commondata
    implicit none
    integer(8) :: k
    real(8) :: f, dfdx
            
    uParticle(k) = -Pi*A*dsin(Pi*f(xParticle(k), timeLocal(k)))*dcos(Pi*yParticle(k))
    vParticle(k) = Pi*A*dcos(Pi*f(xParticle(k), timeLocal(k)))*dsin(Pi*yParticle(k))*dfdx(xParticle(k), timeLocal(k))
        
    xParticle(k+1) = xParticle(k)+uParticle(k)*dt
    yParticle(k+1) = yParticle(k)+vParticle(k)*dt
    
    return
    end subroutine firstEulerUpdate
    
    
    subroutine secondEulerUpdate(k)
    use commondata
    implicit none
    integer(8) :: k
    real(8) :: f, dfdx
    real(8) :: k1x, k1y, k2x, k2y
    
    k1x = -Pi*A*dsin(Pi*f(xParticle(k), timeLocal(k)))*dcos(Pi*yParticle(k))
    k1y =  Pi*A*dcos(Pi*f(xParticle(k), timeLocal(k)))*dsin(Pi*yParticle(k))*dfdx(xParticle(k), timeLocal(k))
    
    k2x = -Pi*A*dsin(Pi*f(xParticle(k)+dt*k1x, timeLocal(k)+dt))*dcos(Pi*(yParticle(k)+dt*k1y))
    k2y =  Pi*A*dcos(Pi*f(xParticle(k)+dt*k1x, timeLocal(k)+dt))*dsin(Pi*(yParticle(k)+dt*k1y))*dfdx(xParticle(k)+dt*k1x, timeLocal(k)+dt)
        
    xParticle(k+1) = xParticle(k)+(k1x+k2x)*dt/2.0d0
    yParticle(k+1) = yParticle(k)+(k1y+k2y)*dt/2.0d0
        
    return
    end subroutine secondEulerUpdate
    
    
    subroutine RK4Update(k)
    use commondata
    implicit none
    integer(8) :: k
    real(8) :: f, dfdx
    real(8) :: k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y
    
    k1x = -Pi*A*dsin(Pi*f(xParticle(k), timeLocal(k)))*dcos(Pi*yParticle(k))
    k1y =  Pi*A*dcos(Pi*f(xParticle(k), timeLocal(k)))*dsin(Pi*yParticle(k))*dfdx(xParticle(k), timeLocal(k))
        
    k2x = -Pi*A*dsin(Pi*f(xParticle(k)+0.5d0*dt*k1x, timeLocal(k)+0.5d0*dt))*dcos(Pi*(yParticle(k)+0.5d0*dt*k1y))
    k2y =  Pi*A*dcos(Pi*f(xParticle(k)+0.5d0*dt*k1x, timeLocal(k)+0.5d0*dt))*dsin(Pi*(yParticle(k)+0.5d0*dt*k1y))*dfdx(xParticle(k)+0.5d0*dt*k1x, timeLocal(k)+0.5d0*dt)
        
    k3x = -Pi*A*dsin(Pi*f(xParticle(k)+0.5d0*dt*k2x, timeLocal(k)+0.5d0*dt))*dcos(Pi*(yParticle(k)+0.5d0*dt*k2y))
    k3y =  Pi*A*dcos(Pi*f(xParticle(k)+0.5d0*dt*k2x, timeLocal(k)+0.5d0*dt))*dsin(Pi*(yParticle(k)+0.5d0*dt*k2y))*dfdx(xParticle(k)+0.5d0*dt*k2x, timeLocal(k)+0.5d0*dt)
        
    k4x = -Pi*A*dsin(Pi*f(xParticle(k)+dt*k3x, timeLocal(k)+dt))*dcos(Pi*(yParticle(k)+dt*k3y))
    k4y =  Pi*A*dcos(Pi*f(xParticle(k)+dt*k3x, timeLocal(k)+dt))*dsin(Pi*(yParticle(k)+dt*k3y))*dfdx(xParticle(k)+dt*k3x, timeLocal(k)+dt)
    
    xParticle(k+1) = xParticle(k)+(k1x+2.0d0*k2x+2.0d0*k3x+k4x)*dt/6.0d0
    yParticle(k+1) = yParticle(k)+(k1y+2.0d0*k2y+2.0d0*k3y+k4y)*dt/6.0d0
        
    return
    end subroutine RK4Update
    
    
    function f(x,t)
    use commondata
    implicit none
    real(8) :: f, x, t
    f = myEpsilon*dsin(omega*t)*x**2.0d0+x-2.0d0*myEpsilon*dsin(omega*t)*x
    
    end function f
    
    
    function dfdx(x,t)
    use commondata
    implicit none
    real(8) :: dfdx, x, t
    dfdx = myEpsilon*dsin(omega*t)*2.0d0*x+1.0d0-2.0d0*myEpsilon*dsin(omega*t)
    
    end function dfdx
    
    
    subroutine output_ASCII()
    use commondata
    implicit none
    integer(8) :: i, j
    character(len=100) :: filename

    write(filename,*) maxTimeStep
    filename = adjustl(filename)

    open(unit=02,file='doubleGyre-'//trim(filename)//'.dat',status='unknown')
    write(02,*) 'TITLE="Lid Driven Cavity"'
    write(02,*) 'VARIABLES="X" "Y" "U" "V" '
    write(02,101) nx, ny
    do j=1,ny
        do i=1,nx
            write(02,100) xGrid(i), yGrid(j), u(i,j), v(i,j)
        enddo
    enddo
100 format(1x,2(e11.4,' '),10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
    close(02)

    return
    end subroutine output_ASCII
    
    
    subroutine LagrangeInterpolation(pointX, pointU, point0X, point0U, order)
    implicit none
    integer(8) :: order
    real(8) :: pointX(1:order+1), pointU(1:order+1)
    real(8) :: point0X, point0U
    real(8) :: TEMP
    integer(8) :: j, k

    point0U = 0.0d0
    do k=1,order+1
        temp = 1.0d0
        do j=1,order+1
            if(j.NE.K) temp=temp*(point0X-pointX(j))/(pointX(k)-pointX(j))
        enddo
        point0U= point0U+temp*pointU(k)
    enddO
                
    return
    end subroutine LagrangeInterpolation
    
    
