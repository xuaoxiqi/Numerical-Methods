#define outputData
!~#define firstEuler
!~#define secondEuler
#define RK4

#define rhsFt
!~#define rhsFty

#define calPosition

    program main
    implicit none
    real(8), parameter :: Pi=4.0d0*datan(1.0d0)
    character*24 ctime, string
    INTEGER*4  time
    real(kind=8) :: start, finish
    integer, parameter :: maxTimeStep=11
    real(8) :: dt
    real(8) :: t(1:maxTimeStep), timeIndex(1:maxTimeStep)
    real(8) :: velocityNumerical(1:maxTimeStep), velocityAnalytical(1:maxTimeStep)
    real(8) :: positionNumerical(1:maxTimeStep), positionAnalytical(1:maxTimeStep)
    real(8) :: acceleration(1:maxTimeStep)
    real(8) :: k1, k2, k3, k4
    integer :: i, interpolationOrder
    integer :: errorNum
    real(8) :: errorL1_1, errorL2_1, errorL1_2, errorL2_2
    
    !~ write(*,*) "Pi=",Pi
    !~ write(*,*) "nint(Pi+0.5)=", nint(Pi+0.5)
    !~ write(*,*) "int(Pi+0.5)=", int(Pi+0.5)
    
    string = ctime( time() )
    write(*,*) 'Start: ', string
    call CPU_TIME(start)
    
#ifdef rhsFt
    write(*,*) "RHS: 2*Pi*cos(2*Pi*t(i))" 
#endif
#ifdef rhsFty
    write(*,*) "RHS: -2*t(i)*u(i)"
#endif

#ifdef calPosition
    write(*,*) "I am calPosition"
#endif
    
    dt = 1.0d0/dble(maxTimeStep-1)
    write(*,*) "dt=", real(dt)
    do i=1,maxTimeStep
        t(i) = (i-1)*dt
        timeIndex(i) = dble(i)
#ifdef rhsFt
        positionAnalytical(i) = -1.0d0/2.0d0/Pi*dcos(2.0d0*Pi*t(i))
        velocityAnalytical(i) = dsin(2.0d0*Pi*t(i))
        acceleration(i) = 2.0d0*Pi*dcos(2.0d0*Pi*t(i))
#endif
#ifdef rhsFty
        velocityAnalytical(i) = dexp(-t(i)**2.0d0)
        acceleration(i) = -2.0d0*t(i)*velocityAnalytical(i)
#endif
    enddo
#ifdef outputData
    open(unit=01,file="analytical.dat",status="unknown")
    do i=1,maxTimeStep
#ifndef calPosition
        write(01,*) t(i), velocityAnalytical(i)
#endif
#ifdef calPosition
        write(01,*) t(i), positionAnalytical(i)
#endif
    enddo
    close(01)
#endif
    
    positionNumerical = 0.0d0
    velocityNumerical = 0.0d0
!---Initial condition
#ifdef rhsFt
    velocityNumerical(1) = 0.0d0
    positionNumerical(1) = -1.0d0/2.0d0/Pi
#endif
#ifdef rhsFty
    velocityNumerical(1) = 1.0d0
#endif
!---Initial condition
    
#ifdef firstEuler
    write(*,*) "I am firstEuler"
    do i=1,maxTimeStep-1
#ifdef rhsFt
        k1 = 2.0d0*Pi*dcos(2.0d0*Pi*t(i))
#endif
#ifdef rhsFty
        k1 = -2.0d0*t(i)*velocityNumerical(i)
#endif
        velocityNumerical(i+1) = velocityNumerical(i)+dt*k1
#ifdef calPosition
        positionNumerical(i+1) = positionNumerical(i)+velocityNumerical(i)*dt  !~!~+0.5d0*dt*dt*k1
#endif
    enddo
#endif
    
#ifdef secondEuler
    write(*,*) "I am secondEuler"
    do i=1,maxTimeStep-1
#ifdef rhsFt
        k1 = 2.0d0*Pi*dcos(2.0d0*Pi*t(i))
        k2 = 2.0d0*Pi*dcos(2.0d0*Pi*t(i+1))
#endif
#ifdef rhsFty
        k1 = -2.0d0*t(i)*velocityNumerical(i)
        k2 = -2.0d0*t(i+1)*(velocityNumerical(i)+dt*k1)
#endif
        velocityNumerical(i+1) = velocityNumerical(i)+dt*(k1+k2)/2.0d0
#ifdef calPosition
        k1 = velocityNumerical(i)
        k2 = velocityNumerical(i+1)
        positionNumerical(i+1) = positionNumerical(i)+dt*(k1+k2)/2.0d0
#endif
    enddo
#endif

#ifdef RK4
    write(*,*) "I am RK4"
!------------------------------------------------------------------
    do i=1,maxTimeStep-1
        interpolationOrder = 5
        if(i.LE.4) interpolationOrder = i
#ifdef rhsFt
        k1 = 2.0d0*Pi*dcos(2.0d0*Pi*t(i))   !~!~2.0d0*Pi*dcos(2.0d0*Pi*t(i)) !~!~t(i)=(i-1)*dt
        k2 = 2.0d0*Pi*dcos(2.0d0*Pi*(t(i)+0.5d0*dt))
        k3 = k2
        k4 = 2.0d0*Pi*dcos(2.0d0*Pi*t(i+1))
#endif
#ifdef rhsFty
        k1 = -2.0d0*t(i)*velocityNumerical(i)
        k2 = -2.0d0*(t(i)+0.5d0*dt)*(velocityNumerical(i)+0.5d0*dt*k1)
        k3 = -2.0d0*(t(i)+0.5d0*dt)*(velocityNumerical(i)+0.5d0*dt*k2)
        k4 = -2.0d0*t(i+1)*(velocityNumerical(i)+dt*k3)
#endif
        velocityNumerical(i+1) = velocityNumerical(i)+dt*(k1+2.0d0*k2+2.0d0*k3+k4)/6.0d0
#ifdef calPosition
        k1 = velocityNumerical(i)
        call LagrangeInterpolation(timeIndex(i-(interpolationOrder-1):i+1), velocityNumerical(i-(interpolationOrder-1):i+1), dble(i)+0.5d0, k2, interpolationOrder)        
        k3 = k2
        k4 = velocityNumerical(i+1)
        positionNumerical(i+1) = positionNumerical(i)+dt*(k1+2.0d0*k2+2.0d0*k3+k4)/6.0d0
#endif
    enddo
#endif

    call CPU_TIME(finish)
    
#ifdef outputData
    open(unit=03,file="numerical.dat",status="unknown")
    do i=1,maxTimeStep
#ifndef calPosition
        write(03,*) t(i), velocityNumerical(i)
#endif
#ifdef calPosition
        write(03,*) t(i), positionNumerical(i)
#endif
    enddo
    close(03)
#endif
    
    errorL1_1 = 0.0d0
    errorL2_1 = 0.0d0
    errorL1_2 = 0.0d0
    errorL2_2 = 0.0d0
    errorNum = 0
    do i=1,maxTimeStep-1
        errorL1_1 = errorL1_1+dabs(velocityNumerical(i+1)-velocityAnalytical(i+1))
        errorL2_1 = errorL2_1+(velocityNumerical(i+1)-velocityAnalytical(i+1))**2.0d0
#ifdef calPosition
        errorL1_2 = errorL1_2+dabs(positionNumerical(i+1)-positionAnalytical(i+1))
        errorL2_2 = errorL2_2+(positionNumerical(i+1)-positionAnalytical(i+1))**2.0d0
#endif
        errorNum = errorNum+1
    enddo
        
    write(*,*) "Time (CPU) = ", real(finish-start), "s"
    write(*,*) "L1 error (velocity)=", errorL1_1/dble(errorNum)
    write(*,*) "L2 error (velocity)=", dsqrt(errorL2_1/dble(errorNum))
#ifdef calPosition
    write(*,*) "L1 error (position)=", errorL1_2/dble(errorNum)
    write(*,*) "L2 error (position)=", dsqrt(errorL2_2/dble(errorNum))
#endif

    string = ctime( time() )
    write(*,*) 'End:   ', string
    
    stop
    end program main
    
    
    subroutine LagrangeInterpolation(pointX, pointU, point0X, point0U, order)
    implicit none
    integer :: order
    real(8) :: pointX(1:order+1), pointU(1:order+1)
    real(8) :: point0X, point0U
    REAL(8) :: TEMP
    integer :: j, k

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