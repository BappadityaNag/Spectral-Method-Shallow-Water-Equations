!----------------------------------------------------------------
! Drivers for Nodal DG 2D Approximation
!----------------------------------------------------------------

program main
    use Constants
    use DGSEM1D_Class
    implicit none

    integer, parameter               :: N = 14   ! Order of Polynomial
    integer, parameter               :: K = 5    ! No. of Elements
    integer, parameter               :: nEqn = 2 ! No. of Physical Variables
    real(kind=rp), dimension(0:K)    :: X
    real(kind=rp), dimension(0:K*N)  :: Y
    type(Mesh1D)                     :: M
    integer                          :: Nt, nstep, i, j
    real(kind=rp)                    :: time, tn, dt
        
    open (unit=10, file='data.txt')

!----------------------------------------------------------------
    
    time = 3.0_rp
    Nt = 1000
    dt = 1.0_rp/Nt
    tn = 0.0_rp

    do i = 0, K
        X(i) = real(i,kind=rp)
    enddo
    M = Construct_Mesh(K,N,nEqn,X)
    write(*,*) "Constructed Mesh ..."
        
    do i = 0, K*N
        Y(i) = real(i,kind=rp)/N
    enddo
    do i = 1, K
        M%e(i)%phi = 0.0_rp
        do j = 0, N
            M%e(i)%phi(j,1) = 2.0_rp*exp( -( ( Y(j+(i-1)*N) - 1.0_rp )/ 0.15_rp )**2.0_rp )
        enddo
    enddo
              
    do nstep = 0, int(Nt*time)
        if (mod(nstep,100) .lt. 1) then
            tn = real(nstep,kind=rp)*dt
            write(*,*) "Timestep : ", tn
            do i = 1, K
                do j = 0, N-1
                    write(10,'(2F20.10)') M%p(i)%eleft-1.0_rp + (M%e(i)%dG%X(j)+1.0_rp)/2.0_rp , M%e(i)%phi(j,1)
                enddo
            enddo
        endif
        call DGSEM1D_by_RK3(tn,dt,M)
    enddo

    call Destruct_Mesh(M)
    write(*,*) "Destroyed Mesh ..."
!----------------------------------------------------------------

end program main
