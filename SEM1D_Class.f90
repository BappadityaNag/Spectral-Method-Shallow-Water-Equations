!----------------------------------------------------------------
! Contains 
! 1. Construct
! 2. Destruct
! 3. Mask
! 4. Unmask
! 5. Global Sum
! 6. Laplace_Operator
! 7. Matrix_Action
!----------------------------------------------------------------

module SEM1D_Class
implicit none
contains

!----------------------------------------------------------------
! Construct the Element
! Usage : gamma = Construct(int N, int K, real {x})
!----------------------------------------------------------------
    function Construct(N, K) result(e)
    
        use Constants
        implicit none

        integer     :: N, K
        type(SEM1D) :: e
        
        e%N = N
        e%K = K
        allocate(e%xi(0:N))
        allocate(e%W(0:N))
        allocate(e%deltax(0:K))
        allocate(e%G(0:N,0:N))
        allocate(e%p(0:K))
        allocate(e%phi(0:N,1:K))
        allocate(e%RHS(0:N,1:K))
        
    end function Construct
!----------------------------------------------------------------

!----------------------------------------------------------------
! Destruct the Element
! Usage : Destruct(SEM1D e)
!----------------------------------------------------------------
    subroutine Destruct(e)

        use Constants
        implicit none

        type(SEM1D)                         :: e

        deallocate(e%xi)
        deallocate(e%W)
        deallocate(e%deltax)
        deallocate(e%G)
        deallocate(e%p)
        deallocate(e%phi)
        deallocate(e%RHS)

    end subroutine Destruct
!----------------------------------------------------------------

!----------------------------------------------------------------
! Global Operations for the 1D-SEM
!----------------------------------------------------------------

! Mask
! Usage : Mask(real {a}, SEM1D e)

    subroutine Mask(a,e)

        use Constants
        implicit none

        real(kind=rp), dimension(0:,1:)     :: a
        type(SEM1D)                         :: e
        integer                             :: k, jR, eR

        do k = 1, e%K-1
            jR = e%p(k)%nodeRight
            eR = e%p(k)%eRight
            a(jR,eR) = 0.0_rp
        enddo

    end subroutine Mask

! UnMask
! Usage : UnMask(real {a}, SEM1D e)
    subroutine UnMask(a,e)
    
        use Constants
        implicit none

        real(kind=rp), dimension(0:,1:)     :: a
        type(SEM1D)                         :: e
        integer                             :: k, jL, jR, eL, eR

        do k = 1, e%K-1
            jR = e%p(k)%nodeRight
            eR = e%p(k)%eRight
            jL = e%p(k)%nodeLeft
            eL = e%p(k)%eLeft
            a(jR,eR) = a(jL,eL)
        enddo

    end subroutine UnMask

! Global_Sum
! Usage : Global_Sum(real {a}, SEM1D e)

    subroutine Global_Sum(a,e)
    
        use Constants
        implicit none

        real(kind=rp), dimension(0:,1:)     :: a
        type(SEM1D)                         :: e
        integer                             :: k, jL, jR, eL, eR
        real(kind=rp)                       :: tmp                                

        do k = 1, e%K-1
            jR = e%p(k)%nodeRight
            eR = e%p(k)%eRight
            jL = e%p(k)%nodeLeft
            eL = e%p(k)%eLeft
            tmp = a(jR,eR) + a(jL,eL)
            a(jR,eR) = tmp
            a(jL,eL) = tmp
        enddo

    end subroutine Global_Sum
!----------------------------------------------------------------

!----------------------------------------------------------------
! Spatial Approximations for the 1D-SEM
!----------------------------------------------------------------

! Laplace_Operator
! Usage : real {D} = Laplace_Operator(real {U}, SEM1D e)

    function Laplace_Operator(U,e) result(D)
    
        use Constants
        use NodalDG_Class
        implicit none

        real(kind=rp), dimension(0:,1:)       :: U
        type(SEM1D)                           :: e
        real(kind=rp), dimension(0:e%N,1:e%K) :: D
        integer                               :: k, j

        do k = 1, e%K
            D(:,k) = MxVDerivative(e%G, U(:,k))
            do j = 0, e%N
                D(j,k) = -2.0_rp*D(j,k)/e%deltax(k)
            enddo
        enddo

    end function Laplace_Operator

! Matrix_Action
! Usage : real {AU} = Matrix_Action(s, real deltat, real {U})

    function Matrix_Action(s,deltat,U,e) result(AU)
    
        use Constants
        implicit none

        real(kind=rp)                         :: s, deltat
        real(kind=rp), dimension(0:,1:)       :: U
        type(SEM1D)                           :: e
        real(kind=rp), dimension(0:e%N,1:e%K) :: AU
        integer                               :: k, j

        call UnMask(U,e)
        AU = Laplace_Operator(U,e)
        do k = 1, e%K
            do j = 0, e%N
                AU(j,k) = e%W(j)/2.0_rp*e%deltax(k)*U(j,k) + s*deltat/2.0_rp*AU(j,k)
            enddo
        enddo
        call Global_Sum(AU,e)
        call Mask(U,e)          ! Why is this line required?
        if (s .lt. 0.0_rp) then ! for the LHS
            call Mask(AU,e)
            AU(0,1) = 0.0_rp
            AU(e%N,e%K) = 0.0_rp
        endif

    end function Matrix_Action

! Residual
! Usage : real {r} = Residual(SEM1D e)

    function Residual(e, deltat) result(r)
    
        use Constants
        implicit none

        type(SEM1D)                           :: e
        real(kind=rp), dimension(0:e%N,1:e%K) :: r
        real(kind=rp)                         :: deltat
        integer                               :: k, j

        r = Matrix_Action(1.0_rp,deltat,e%phi,e)
        do k = 1, e%K
            do j = 0, e%N
                r(j,k) = e%RHS(j,k) - r(j,k)
            enddo
        enddo
        call Mask(r,e)
        r(0,1) = 0.0_rp
        r(e%N,e%K) = 0.0_rp

    end function Residual
!----------------------------------------------------------------

end module SEM1D_Class
