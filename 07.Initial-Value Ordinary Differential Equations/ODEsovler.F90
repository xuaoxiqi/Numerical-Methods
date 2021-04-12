#define outputData
!~#define firstEuler
!~#define secondEuler
#define RK4

    program main
    implicit none
    real(8), parameter :: Pi=4.0d0*datan(1.0d0)
    character*24 ctime, string
    INTEGER*4  time
    real(kind=8) :: start, finish
    integer, parameter :: meshMax=11
    real(8) :: dt
    real(8) :: t(1:meshMax), tMesh(1:meshMax)
    real(8) :: velocityNumerical(1:meshMax), velocityAnalytical(1:meshMax)
    real(8) :: acceleration(1:meshMax)
    real(8) :: k1, k2, k3, k4
    integer :: i
    integer :: errorNum
    real(8) :: errorL1, errorL2
    
    write(*,*) "Pi=",Pi
    write(*,*) "nint(Pi+0.5)=", nint(Pi+0.5)
    write(*,*) "int(Pi+0.5)=", int(Pi+0.5)
    
    string = ctime( time() )
    write(*,*) 'Start: ', string
    call CPU_TIME(start)
    
    dt = 1.0/dble(meshMax-1)
    do i=1,meshMax
        t(i) = (i-1)*dt
        acceleration(i) = 2.0d0*Pi*dcos(2.0d0*Pi*t(i))
        velocityAnalytical(i) = dsin(2.0d0*Pi*t(i))
        tMesh(i) = dble(i)
    enddo
#ifdef outputData
    open(unit=01,file="analytical.dat",status="unknown")
    do i=1,meshMax
        write(01,*) t(i), velocityAnalytical(i)
    enddo
    close(01)
#endif
    
    velocityNumerical = 0.0d0
    errorL1 = 0.0d0
    errorL2 = 0.0d0
    errorNum = 0
#ifdef firstEuler    
    do i=1,meshMax-1
        velocityNumerical(i+1) = velocityNumerical(i)+dt*acceleration(i)
        
        errorL1 = errorL1+dabs(velocityNumerical(i+1)-velocityAnalytical(i+1))
        errorL2 = errorL2+(velocityNumerical(i+1)-velocityAnalytical(i+1))**2.0d0
        errorNum = errorNum+1
    enddo
#endif
    
#ifdef secondEuler
    do i=1,meshMax-1
        velocityNumerical(i+1) = velocityNumerical(i)+dt*(acceleration(i)+acceleration(i+1))/2.0d0

        errorL1 = errorL1+dabs(velocityNumerical(i+1)-velocityAnalytical(i+1))
        errorL2 = errorL2+(velocityNumerical(i+1)-velocityAnalytical(i+1))**2.0d0
        errorNum = errorNum+1
    enddo
#endif

#ifdef RK4
    do i=1,3
        velocityNumerical(i+1) = velocityNumerical(i)+dt*(acceleration(i)+acceleration(i+1))/2.0d0
        
        errorL1 = errorL1+dabs(velocityNumerical(i+1)-velocityAnalytical(i+1))
        errorL2 = errorL2+(velocityNumerical(i+1)-velocityAnalytical(i+1))**2.0d0
        errorNum = errorNum+1
    enddo
    
    do i=4,meshMax-1
        k1 = acceleration(i)
        call LagrangeInterpolation(tMesh(i-3:i+1), acceleration(i-3:i+1), dble(i)+0.5d0, k2, 3)
        k3 = k2
        k4 = acceleration(i+1)
        velocityNumerical(i+1) = velocityNumerical(i)+dt*(k1+2.0d0*k2+2.0d0*k3+k4)/6.0d0
        
        errorL1 = errorL1+dabs(velocityNumerical(i+1)-velocityAnalytical(i+1))
        errorL2 = errorL2+(velocityNumerical(i+1)-velocityAnalytical(i+1))**2.0d0
        errorNum = errorNum+1
    enddo
#endif

#ifdef outputData
    open(unit=03,file="numerical.dat",status="unknown")
    do i=1,meshMax
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