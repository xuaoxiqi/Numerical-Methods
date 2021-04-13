#define outputData
!~#define firstEuler
!~#define secondEuler
#define RK4

!~#define rhsFt
#define rhsFty

    program main
    implicit none
    real(8), parameter :: Pi=4.0d0*datan(1.0d0)
    character*24 ctime, string
    INTEGER*4  time
    real(kind=8) :: start, finish
    integer, parameter :: maxTimeStep=11
    real(8) :: dt
    real(8) :: t(1:maxTimeStep)
    real(8) :: velocityNumerical(1:maxTimeStep), velocityAnalytical(1:maxTimeStep)
    real(8) :: acceleration(1:maxTimeStep)
    real(8) :: k1, k2, k3, k4
    integer :: i
    integer :: errorNum
    real(8) :: errorL1, errorL2
    
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
    
    dt = 1.0d0/dble(maxTimeStep-1)
    write(*,*) "dt=", dt
    do i=1,maxTimeStep
        t(i) = (i-1)*dt
#ifdef rhsFt
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
        write(01,*) t(i), velocityAnalytical(i)
    enddo
    close(01)
#endif
    
    velocityNumerical = 0.0d0
    errorL1 = 0.0d0
    errorL2 = 0.0d0
    errorNum = 0
!---Initial condition
#ifdef rhsFt
    velocityNumerical(1) = 0.0d0
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
        
        errorL1 = errorL1+dabs(velocityNumerical(i+1)-velocityAnalytical(i+1))
        errorL2 = errorL2+(velocityNumerical(i+1)-velocityAnalytical(i+1))**2.0d0
        errorNum = errorNum+1
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

        errorL1 = errorL1+dabs(velocityNumerical(i+1)-velocityAnalytical(i+1))
        errorL2 = errorL2+(velocityNumerical(i+1)-velocityAnalytical(i+1))**2.0d0
        errorNum = errorNum+1
    enddo
#endif

#ifdef RK4
    write(*,*) "I am RK4"
    do i=1,maxTimeStep-1
    
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
        
        errorL1 = errorL1+dabs(velocityNumerical(i+1)-velocityAnalytical(i+1))
        errorL2 = errorL2+(velocityNumerical(i+1)-velocityAnalytical(i+1))**2.0d0
        errorNum = errorNum+1
    enddo
#endif

#ifdef outputData
    open(unit=03,file="numerical.dat",status="unknown")
    do i=1,maxTimeStep
        write(03,*) t(i), velocityNumerical(i)
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
    
    