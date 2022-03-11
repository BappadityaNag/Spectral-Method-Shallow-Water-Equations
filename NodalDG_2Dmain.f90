!----------------------------------------------------------------
! Drivers for Nodal DG 2D Approximation
!----------------------------------------------------------------

program main
    use Constants
    use NodalDG_2D_Class
    implicit none

    type(NodalDG_2DSystem)               :: gamma
    integer, parameter                   :: N = 10, M = 5, nEqn = 3
    real(kind=rp)                        :: c = 10.0_rp    ! Phase Speed
    integer                              :: Nt, nstep, i, j
    real(kind=rp), dimension(0:N)        :: X
    real(kind=rp), dimension(0:M)        :: Y
    real(kind=rp)                        :: time, tn, dt
        
    open (unit=10, file='data.txt')

!----------------------------------------------------------------
    
    gamma = Construct(nEqn,N,M)
        
    time = 0.5_rp
    Nt = 100
    dt = time/Nt
    tn = 0.0_rp

    do i = 0, N
        X(i) = -1.0_rp + 2.0_rp*real(i,kind=rp)/N
    enddo
    do j = 0, M
        Y(j) = -1.0_rp + 2.0_rp*real(j,kind=rp)/M
    enddo
    gamma%Q = 0.0_rp
    
    do i = 0, N
        do j = 0, M
            ! Circular Wave
            gamma%Q(i,j,1) = 0.5*exp(-10.*(X(i)**2+Y(j)**2))
            ! Benchmark Circular Wave
            ! gamma%Q(i,j,1) = exp(-log(2.0_rp)/(0.06)**2*(X(i)**2+Y(j)**2))
        enddo
    enddo

    do nstep = 0, Nt-1
        write(*,*) "Timestep : ", nstep+1
        call DG2DStep_by_RK3(tn,dt,gamma,c)
        tn = (real(nstep,kind=rp)+1.0_rp)*dt
        do i = 0, N
            write(10,'(6F20.15)') gamma%Q(i,:,1)
        enddo
    enddo

    call Destruct(gamma)
!----------------------------------------------------------------

end program main
