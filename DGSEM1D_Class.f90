!-----------------------------------------------------------------------
! Contains 
!  1. Construct_Element
!  2. Destruct_Element
!  3. Construct_Mesh
!  4. Destruct_Mesh
!  5. DGSEM1D_by_RK3
!  6. Global_Time_Derivative
!  7. Riemann_Solver
!  8. Flux_DGSEM1D
!  9. External_State
! 10. Local_Time_Derivative
! 11. Interpolate_to_Boundaries
! 12. Affine_Map
!-----------------------------------------------------------------------

module DGSEM1D_Class
implicit none
contains

!-----------------------------------------------------------------------
! Construct the Element
! Usage : E = Construct_Element(N, nEqn, dG, xL, xR)
!-----------------------------------------------------------------------
    function Construct_Element(N,nEqn,dG,xL,xR) result(E)

        use Constants
        implicit none

        type(NodalDG_Storage)          :: dG
        integer                        :: N      ! No. of nodal points 
        integer                        :: nEqn   ! No. of eqns in the physical system
        real(kind=rp)                  :: xL, xR ! Size and left and right boundaries of the element
        type(Element)                  :: E

        E%N = N
        E%nEqn = nEqn
        E%xL = xL
        E%xR = xR
        E%deltax = xR-xL
        E%dG = dG
        allocate(E%phi(0:N,1:nEqn))    ! Solution
        allocate(E%phidot(0:N,1:nEqn)) ! Time derivative of the Solution
        allocate(E%G(0:N,1:nEqn))      ! For low-storage Runga-Kutta
        allocate(E%phiL(1:nEqn))       ! Solution on left element boundary
        allocate(E%phiR(1:nEqn))       ! Solution on right element boundary
        allocate(E%FL(1:nEqn))         ! Flux on left element boundary
        allocate(E%FR(1:nEqn))         ! Flux on right element boundary
        
    end function Construct_Element
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Destruct the Element
! Usage : Destruct_Element(E)
!-----------------------------------------------------------------------
    subroutine Destruct_Element(E)

        use Constants
        implicit none

        type(Element)                  :: E

        deallocate(E%phi)              ! Solution
        deallocate(E%phidot)           ! Time derivative of the Solution
        deallocate(E%G)                ! For low-storage Runga-Kutta
        deallocate(E%phiL)             ! Solution on left element boundary
        deallocate(E%phiR)             ! Solution on right element boundary
        deallocate(E%FL)               ! Flux on left element boundary
        deallocate(E%FR)               ! Flux on right element boundary
        
    end subroutine Destruct_Element
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Construct_Mesh
! Usage : M = Construct_Mesh(K, N, real {x})
!-----------------------------------------------------------------------
    function Construct_Mesh(K, N, nEqn, x) result(M)
    
        use Constants
        use NodalDG_Class
        implicit none

        integer                        :: K, N, nEqn, j
        real(kind=rp), dimension(0:)   :: x
        type(NodalDG_Storage)          :: dG
        type(Mesh1D)                   :: M

        dG = Construct(N)
        M%K = K
        allocate(M%e(1:K))
        do j = 1, M%K
            M%e(j) = Construct_Element(N, nEqn, dG, x(j-1), x(j))
        enddo
        call Destruct(dG)
        
        allocate(M%p(0:M%K))
        do j = 1, M%K-1
            M%p(j)%eLeft  = j
            M%p(j)%eRight = j + 1
        enddo
        M%p(0)%eLeft  = Null
        M%p(0)%eRight = 1
        M%p(K)%eLeft  = K
        M%p(K)%eRight = Null

    end function Construct_Mesh
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Usage : call Destruct_Mesh(M)
!-----------------------------------------------------------------------
    subroutine Destruct_Mesh(M)
    
        use Constants
        implicit none

        type(Mesh1D)                   :: M
        integer                        :: k

        do k = 1, M%K
             call Destruct_Element(M%e(k))
        enddo
        deallocate(M%e)
        deallocate(M%p)
        
    end subroutine Destruct_Mesh
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Global_Mesh_Procedures : Mesh Global Procedures for the DGSEM
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Low Storage Runge-Kutta Integration of DGSEM 1D Approximation
! Usage : DGSEM1D_by_RK3(P,D,tn,dt)
!-----------------------------------------------------------------------
    subroutine DGSEM1D_by_RK3(tn,dt,M)
        use Constants
        implicit none

        type(Mesh1D)                   :: M
        real(kind=rp)                  :: t, tn, dt
        integer                        :: i, j, k, l

        real(kind=rp), dimension(3) :: a = (/ 0.0_rp, -5.0_rp/9.0_rp, -153.0_rp/128.0_rp /)
        real(kind=rp), dimension(3) :: b = (/ 0.0_rp, 1.0_rp/3.0_rp, 3.0_rp/4.0_rp /)
        real(kind=rp), dimension(3) :: g = (/ 1.0_rp/3.0_rp, 15.0_rp/16.0_rp, 8.0_rp/15.0_rp /)

        do k = 1, 3
            t = tn + b(k)*dt
            call Global_Time_Derivative(M)
            do i = 1, M%K
                do j = 0, M%e(i)%N
                    do l = 1, M%e(i)%nEqn
                        M%e(i)%G(j,l) = a(k)*M%e(i)%G(j,l) + M%e(i)%phiDot(j,l)
                        M%e(i)%phi(j,l) = M%e(i)%phi(j,l) + g(k)*dt*M%e(i)%G(j,l)
                    enddo
                enddo
            enddo
        enddo

    end subroutine DGSEM1D_by_RK3
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Global Time Derivative
! Usage : Global_Time_Derivative(t,M,c)
!-----------------------------------------------------------------------
    subroutine Global_Time_Derivative(M)
    
        use Constants
        use NodalDG_2D_Class
        implicit none

        type(Mesh1D)                               :: M
        integer                                    :: k, idL, idR
        real(kind=rp), dimension(1:M%e(1)%nEqn)    :: QextL, QextR
        real(kind=rp), dimension(1:M%e(1)%nEqn)    :: FL, FR         ! Flux on left/right element boundary
        
        do k = 1, M%K
            call Interpolate_to_Boundaries(M%e(k))
        enddo
        
        k  = M%p(0)%eRight
        QextL = External_State( M%e(k)%phiL, Affine_Map( M%e(k)%xL,M%e(k)%deltax,-1.0_rp ) )
        k  = M%p(M%K)%eLeft
        QextR = External_State( M%e(k)%phiR, Affine_Map( M%e(k)%xL,M%e(k)%deltax,+1.0_rp ) )
        
        do k = 0, M%K
            idL = M%p(k)%eLeft
            idR = M%p(k)%eRight
            if (idL .eq. Null) then
                M%e(idR)%FL = Riemann_Solver( QextL, M%e(idR)%phiL, -1.0_rp )
            elseif (idR .eq. Null) then
                M%e(idL)%FR = Riemann_Solver( M%e(idL)%phiR, QextR, +1.0_rp )
            else
                FL = Riemann_Solver( M%e(idL)%phiR, M%e(idR)%phiL, +1.0_rp )
                M%e(idR)%FL = - FL
                M%e(idL)%FR =   FL
            endif
        enddo
        do k = 1, M%K
            call Local_Time_Derivative(M%e(k))
        enddo

    end subroutine Global_Time_Derivative
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Riemann Solver
! Usage : real {Fstar} = Riemann_Solver(real {QL}, real {QR}, normal)
!-----------------------------------------------------------------------

    function Riemann_Solver(QL, QR, normal) result(Fstar)
        use Constants
        implicit none

        real(kind=rp), dimension(1:)             :: QL, QR
        real(kind=rp)                            :: normal
        real(kind=rp), dimension(1:ubound(QL,1)) :: Fstar
        real(kind=rp)                :: uL, uR, vL, vR, wplus, wminus

         uL = QL(1); vL = QL(2)
         uR = QR(1); vR = QR(2)
         wplus    = uL+vL
         wminus   = uR-vR
         Fstar(1) = 0.5_rp*(wplus-wminus)*normal
         Fstar(2) = 0.5_rp*(wplus+wminus)*normal
         !Fstar(1) = 0.5_rp*(3.0_rp*wplus-wminus)*normal
         !Fstar(2) = 0.5_rp*(3.0_rp*wplus+wminus)*normal

    end function Riemann_Solver
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Flux Vector for the SEM 1D Wave Equation
! Usage : real {F} = Flux_DGSEM1D(real {Q}, real c)
!-----------------------------------------------------------------------
    function Flux_DGSEM1D(Q) result(F)
        use Constants
        implicit none

        real(kind=rp), dimension(1:)            :: Q
        real(kind=rp), dimension(1:ubound(Q,1)) :: F

        F(1) = Q(2)
        F(2) = Q(1)
        !F(1) = Q(1) + 2.0_rp*Q(2)
        !F(2) = 2.0_rp*Q(1) + Q(2)
        
    end function Flux_DGSEM1D
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Boundary Conditions
! Usage : {Qex} = External_State( real {Qin}, -1, t, RIGHT )
!-----------------------------------------------------------------------
    function External_State( Qin, x ) result(Qex)
        use Constants
        implicit none

        real(kind=rp), dimension(1:)              :: Qin
        real(kind=rp), dimension(1:ubound(Qin,1)) :: Qex
        real(kind=rp)                             :: x
        
        if (x .eq. 0.0_rp) then 
            Qex(1) =  Qin(1)
            Qex(2) = -Qin(2)
        elseif (x .eq. 5.0_rp) then
            Qex = Qin
            !Qex = 0.0_rp
        endif
        
    end function External_State
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Local_DSEM_Procedures : Local Procedures for the DGSEM
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Local Time Derivative
! Usage : call Local_Time_Derivative(Element E)
!-----------------------------------------------------------------------
    subroutine Local_Time_Derivative(E)

        use Constants
        use NodalDG_2D_Class
        implicit none

        integer                                  :: j, k
        type(Element)                            :: E
        real(kind=rp), dimension(0:E%N,1:E%nEqn) :: F, Fderiv

        do j = 0, E%dG%N
            F(j,:) = Flux_DGSEM1D(E%phi(j,:))
        enddo
        Fderiv = SystemDG_Derivative(E%FL,E%FR,F,E%dG%D1,E%dG%lleft,E%dG%lright,E%dG%W)
        do j = 0, E%dG%N
            do k = 1, E%nEqn
                E%phidot(j,k) = -2.0_rp*Fderiv(j,k)/E%deltax
            enddo
        enddo
        
    end subroutine Local_Time_Derivative
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Interpolate to Boundaries
! Usage : call Interpolate_to_Boundaries(Element E)
!-----------------------------------------------------------------------
    subroutine Interpolate_to_Boundaries(E)

        use Constants
        use NodalDG_2D_Class
        implicit none

        integer                        :: k
        type(Element)                  :: E

        do k = 1, E%nEqn
            E%phiL(k) = Interpolate_to_Boundary(E%phi(:,k),E%dG%lleft)
            E%phiR(k) = Interpolate_to_Boundary(E%phi(:,k),E%dG%lright)
        enddo
        
    end subroutine Interpolate_to_Boundaries
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Affine Map
! Usage : real x = Affine_Map(xleft, deltax, xi)
!-----------------------------------------------------------------------
    function Affine_Map(xleft, deltax, xi) result(x)

        use Constants
        implicit none

        real(kind=rp)          :: xleft, deltax, xi, x

        x = xleft + (xi+1.0_rp)/2*deltax
        
    end function Affine_Map
!-----------------------------------------------------------------------

end module DGSEM1D_Class
