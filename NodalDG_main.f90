!----------------------------------------------------------------
! Drivers for Nodal DG Approximation
!----------------------------------------------------------------

program main
    use Constants
    use NodalDG_Class
    implicit none

    type(NodalDG_Storage)                :: gamma
    integer, parameter                   :: N = 20
    real(kind=rp)                        :: c =  -30.0_rp     ! Phase Speed
    real(kind=rp)                        :: bc = 0.0_rp    ! Boundary Condition
    integer                              :: Nt, nstep, j
    real(kind=rp)                        :: time, tn, dt
        
    open (unit=10, file='data.txt')

!----------------------------------------------------------------
    
    gamma = Construct(N)
        
    time = 0.1_rp
    Nt = 1000
    dt = time/Nt
    tn = 0.0_rp
    
    do j = 0, N
        gamma%X(j) = -1.0_rp + 2.0_rp*real(j,kind=rp)/N
        gamma%phi(j) = 0.5*exp(-10.*(gamma%X(j))**2)
        !gamma%phi(j) = 0.5
    enddo

    do nstep = 0, Nt-1
        call DGStep_by_RK3(tn,dt,gamma,c,bc)
        tn = (real(nstep,kind=rp)+1.0_rp)*dt
        write(10,'(21F20.15)') gamma%phi
    enddo

    call Destruct(gamma)

!----------------------------------------------------------------

end program main
