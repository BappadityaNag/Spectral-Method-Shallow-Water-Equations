!----------------------------------------------------------------
! Contains 
! 1. Construct
! 2. Destruct
! 3. DG2DStep_by_RK3
! 4. SystemDG_Derivative
! 5. DG2D_Time_Derivative
! 6. Riemann_Solver
! 7. Wave_Equation_Fluxes (xFlux, yFlux)
! 8. MxVDerivative
!----------------------------------------------------------------

module NodalDG_2D_Class
implicit none
contains

!----------------------------------------------------------------
! Construct the Contour
! Usage : gamma = Construct(nEqn, N, M)
!----------------------------------------------------------------
    function Construct(nEqn, N, M) result(gamma)
        use Constants
        use Quadrature_Points_and_Weights_Class
        use Lagrange_Interpolation_Class
        implicit none

        integer                                    :: N, M, nEqn
        type(NodalDG_2DSystem)                     :: gamma
        integer                                    :: i, j
        real(kind=rp), dimension(:,:), allocatable :: D

        gamma%nEqn = nEqn
        gamma%spA%N = N
        allocate(gamma%Q(0:N,0:M,1:nEqn))
        
        allocate(gamma%spA%xi(0:N))
        allocate(gamma%spA%Wxi(0:N))
        allocate(gamma%spA%WBxi(0:N))
        allocate(gamma%spA%D1xi(0:N,0:N))
        allocate(gamma%spA%D2xi(0:N,0:N))
        allocate(gamma%spA%lxileft(0:N))
        allocate(gamma%spA%lxiright(0:N))
        allocate(D(0:N,0:N))
        
        call Legendre_Gauss_Nodes_and_Weights( N,gamma%spA%xi,gamma%spA%Wxi )
        gamma%spA%WBxi = Barycentric_Weights( gamma%spA%xi )
        gamma%spA%lxileft = Lagrange_Interpolating_Polynomials( -1.0_rp,gamma%spA%xi,gamma%spA%WBxi )
        gamma%spA%lxiright = Lagrange_Interpolating_Polynomials( 1.0_rp,gamma%spA%xi,gamma%spA%WBxi )
        D = Polynomial_Derivative_Matrix(gamma%spA%xi)
        
        do j = 0,N
            do i = 0,N
                gamma%spA%D1xi(i,j) = -D(j,i)*gamma%spA%Wxi(j)/gamma%spA%Wxi(i)
            enddo
        enddo
        
        deallocate(D)

        gamma%spA%M = M
        allocate(gamma%spA%eta(0:M))
        allocate(gamma%spA%Weta(0:M))
        allocate(gamma%spA%WBeta(0:M))
        allocate(gamma%spA%D1eta(0:M,0:M))
        allocate(gamma%spA%D2eta(0:M,0:M))
        allocate(gamma%spA%letaleft(0:M))
        allocate(gamma%spA%letaright(0:M))
        allocate(D(0:M,0:M))

        call Legendre_Gauss_Nodes_and_Weights( M,gamma%spA%eta,gamma%spA%Weta )
        gamma%spA%WBeta = Barycentric_Weights( gamma%spA%eta )
        gamma%spA%letaleft = Lagrange_Interpolating_Polynomials( -1.0_rp,gamma%spA%eta,gamma%spA%WBeta )
        gamma%spA%letaright = Lagrange_Interpolating_Polynomials( 1.0_rp,gamma%spA%eta,gamma%spA%WBeta )
        D = Polynomial_Derivative_Matrix(gamma%spA%eta)
        
        do j = 0,M
            do i = 0,M
                gamma%spA%D1eta(i,j) = -D(j,i)*gamma%spA%Weta(j)/gamma%spA%Weta(i)
            enddo
        enddo

        deallocate(D)

    end function Construct
!----------------------------------------------------------------

!----------------------------------------------------------------
! Destruct the Contour
! Usage : Destruct({contour})
!----------------------------------------------------------------
    subroutine Destruct(gamma)
        use Constants
        implicit none

        type(NodalDG_2DSystem) :: gamma

        deallocate(gamma%Q)

        deallocate(gamma%spA%xi)
        deallocate(gamma%spA%Wxi)
        deallocate(gamma%spA%WBxi)
        deallocate(gamma%spA%D1xi)
        deallocate(gamma%spA%D2xi)
        deallocate(gamma%spA%lxileft)
        deallocate(gamma%spA%lxiright)

        deallocate(gamma%spA%eta)
        deallocate(gamma%spA%Weta)
        deallocate(gamma%spA%WBeta)
        deallocate(gamma%spA%D1eta)
        deallocate(gamma%spA%D2eta)
        deallocate(gamma%spA%letaleft)
        deallocate(gamma%spA%letaright)

    end subroutine Destruct
!----------------------------------------------------------------

!----------------------------------------------------------------
! Low Storage Runge-Kutta Integration of a Nodal DG 2D Approximation
! Usage : DG2DStep_by_RK3(P,D,tn,dt)
!----------------------------------------------------------------
    subroutine DG2DStep_by_RK3(tn,dt,gamma,c)
        use Constants
        implicit none

        type(NodalDG_2DSystem)                 :: gamma
        real(kind=rp), dimension(0:ubound(gamma%Q,1),0:ubound(gamma%Q,2),1:gamma%nEqn) :: Qdot, tmp
        real(kind=rp)                          :: t, tn, dt, c
        integer                                :: N, M, nEqn, i, j, k, l

        real(kind=rp), dimension(3) :: a = (/ 0.0_rp, -5.0_rp/9.0_rp, -153.0_rp/128.0_rp /)
        real(kind=rp), dimension(3) :: b = (/ 0.0_rp, 1.0_rp/3.0_rp, 3.0_rp/4.0_rp /)
        real(kind=rp), dimension(3) :: g = (/ 1.0_rp/3.0_rp, 15.0_rp/16.0_rp, 8.0_rp/15.0_rp /)

        N    = ubound(gamma%Q,1)
        M    = ubound(gamma%Q,2)
        nEqn = gamma%nEqn
        tmp = 0.0_rp
        do k = 1,3
            t = tn + b(k)*dt
            QDot = DG2D_Time_Derivative(t,gamma,c)
            do i = 0,N
                do j = 0,M
                    do l = 1,nEqn
                        tmp(i,j,l) = a(k)*tmp(i,j,l) + QDot(i,j,l)
                        gamma%Q(i,j,l) = gamma%Q(i,j,l) + g(k)*dt*tmp(i,j,l)
                    enddo
                enddo
            enddo
        enddo

    end subroutine DG2DStep_by_RK3
!----------------------------------------------------------------

!----------------------------------------------------------------
! Time Derivative in 2D using the Discontinuous Galerkin Approximation
! Usage : {Qdot} = DG2D_Time_Derivative(t)
!----------------------------------------------------------------
    function DG2D_Time_Derivative(t,gamma,c) result(Qdot)
        use Constants
        implicit none

        type(NodalDG_2DSystem)                 :: gamma
        real(kind=rp)                          :: t, c, x, y
        real(kind=rp), dimension(1:gamma%nEqn) :: QintL, QintR, QextL, QextR
        real(kind=rp), dimension(1:gamma%nEqn) :: FstarL, FstarR, GstarL, GstarR
        real(kind=rp), dimension(0:ubound(gamma%Q,1),1:gamma%nEqn) :: F, Fderiv
        real(kind=rp), dimension(0:ubound(gamma%Q,2),1:gamma%nEqn) :: G, Gderiv
        real(kind=rp), dimension(0:ubound(gamma%Q,1),0:ubound(gamma%Q,2),1:gamma%nEqn) :: Qdot
        integer                                :: N, M, nEqn, i, j, k
        
        N    = gamma%spA%N
        M    = gamma%spA%M
        nEqn = gamma%nEqn

        do j = 0, M
            y = gamma%spA%eta(j)
            do k = 1, nEqn
                QintL(k) = Interpolate_to_Boundary( gamma%Q(:,j,k), gamma%spA%lxileft )
                QintR(k) = Interpolate_to_Boundary( gamma%Q(:,j,k), gamma%spA%lxiright )
            enddo
            QextL = External_State( QintL, -1.0_rp, y ) ! t, LEFT
            QextR = External_State( QintR, 1.0_rp, y )  ! t, RIGHT
            FstarL = Riemann_Solver( QintL, QextL, (/-1.0_rp,0.0_rp/), c )
            FstarR = Riemann_Solver( QintR, QextR, (/1.0_rp,0.0_rp/), c )
            do i = 0, N
                F(i,:) = xFlux( gamma%Q(i,j,:), c )
            enddo
            Fderiv = SystemDG_Derivative( FstarL,FstarR,F,gamma%spA%D1xi,gamma%spA%lxileft,gamma%spA%lxiright,gamma%spA%Wxi )
            do i = 0, N
                do k = 1, nEqn
                    Qdot(i,j,k) = - Fderiv(i,k)
                enddo
            enddo
        enddo

        do i = 0, N
            x = gamma%spA%xi(j)
            do k = 1, nEqn
                QintL(k) = Interpolate_to_Boundary( gamma%Q(i,:,k), gamma%spA%letaleft )
                QintR(k) = Interpolate_to_Boundary( gamma%Q(i,:,k), gamma%spA%letaright )
            enddo
            QextL = External_State( QintL, x, -1.0_rp ) ! t, BOTTOM
            QextR = External_State( QintR, x, 1.0_rp )  ! t, TOP
            GstarL = Riemann_Solver( QintL, QextL, (/0.0_rp,-1.0_rp/), c )
            GstarR = Riemann_Solver( QintR, QextR, (/0.0_rp,1.0_rp/), c )
            do j = 0, M
                G(j,:) = yFlux( gamma%Q(i,j,:), c )
            enddo
            Gderiv = SystemDG_Derivative( GstarL,GstarR,G,gamma%spA%D1eta,gamma%spA%letaleft,gamma%spA%letaright,gamma%spA%Weta )
            do j = 0, M
                do k = 1, nEqn
                    Qdot(i,j,k) = Qdot(i,j,k) - Gderiv(j,k)
                enddo
            enddo
        enddo

    end function DG2D_Time_Derivative
!----------------------------------------------------------------

!----------------------------------------------------------------
! First Derivative using the Discontinuous Galerkin Approximation
! Usage : {Fderiv} = SystemDG_Derivative({FL},{FR},{F},{D},{lleft},{lright},{W})
!----------------------------------------------------------------
    function SystemDG_Derivative(FL,FR,F,D,lleft,lright,W) result(Fderiv)
        use Constants
        implicit none

        real(kind=rp), dimension(1:)                          :: FL, FR
        real(kind=rp), dimension(0:)                          :: lleft, lright, W
        real(kind=rp), dimension(0:,0:)                       :: D
        real(kind=rp), dimension(0:,1:)                       :: F
        real(kind=rp), dimension(0:ubound(F,1),1:ubound(F,2)) :: Fderiv
        integer                                               :: N, nEqn, j, k

        N    = ubound(F,1)
        nEqn = ubound(F,2)
        
        do k = 1, nEqn
            Fderiv(:,k) = MxVDerivative(D,F(:,k))
        enddo
        do j = 0, N
            do k = 1, nEqn
                Fderiv(j,k) = Fderiv(j,k) + ( FR(k)*lright(j) + FL(k)*lleft(j) ) / W(j)
            enddo
        enddo
        
    end function SystemDG_Derivative
!----------------------------------------------------------------

!----------------------------------------------------------------
! Boundary Conditions
! Usage : {Qex} = External_State( real {Qin}, -1, y, t, LEFT )
!----------------------------------------------------------------
    function External_State( Qin, x, y ) result(Qex)
        use Constants
        implicit none

        real(kind=rp), dimension(1:)              :: Qin
        real(kind=rp), dimension(1:ubound(Qin,1)) :: Qex
        real(kind=rp)                             :: x, y
        
        if (x .eq. 1.0_rp) then 
            Qex(1) =   Qin(1)
            Qex(2) = - Qin(2)
            Qex(3) =   Qin(3)
        else
            Qex = 0.0_rp
        endif
        
    end function External_State
!----------------------------------------------------------------

!----------------------------------------------------------------
! Wave_Equation_Fluxes
! Flux Vectors for the Two Dimensional Wave Equation
! Usage : real {F} = xFlux(real {Q}, real c)
!----------------------------------------------------------------

    function xFlux(Q,c) result(F)
        use Constants
        implicit none

        real(kind=rp), dimension(1:)                          :: Q
        real(kind=rp), dimension(1:ubound(Q,1))               :: F
        real(kind=rp)                                         :: c

        F(1) = c**2 * Q(2)
        F(2) = Q(1)
        F(3) = 0.0_rp
    end function xFlux

    function yFlux(Q,c) result(G)
        use Constants
        implicit none

        real(kind=rp), dimension(1:)                          :: Q
        real(kind=rp), dimension(1:ubound(Q,1))               :: G
        real(kind=rp)                                         :: c

        G(1) = c**2 * Q(3)
        G(2) = 0.0_rp
        G(3) = Q(1)
    end function yFlux
    
!----------------------------------------------------------------

!!----------------------------------------------------------------
!! Riemann Solver
!! Usage : real {Fstar} = Riemann_Solver(real {QL}, real {QR}, normal)
!!----------------------------------------------------------------

!    function Riemann_Solver(QL, QR, normal, c) result(Fstar)
!        use Constants
!        implicit none

!        real(kind=rp), dimension(1:)             :: QL, QR
!        real(kind=rp), dimension(1:2)            :: normal
!        real(kind=rp), dimension(1:ubound(QL,1)) :: Fstar
!        real(kind=rp)        :: c, pL, pR, uL, uR, vL, vR, wplus, wminus

!        pL = QL(1); uL = QL(2); vL = QL(3)
!        pR = QR(1); uR = QR(2); vR = QR(3)
!        wplus    = pL + c*( normal(1)*uL + normal(2)*vL )
!        wminus   = pR - c*( normal(1)*uR + normal(2)*vR )
!        Fstar(1) = c*(wplus-wminus)/2.0_rp
!        Fstar(2) = normal(1)*(wplus+wminus)/2.0_rp
!        Fstar(3) = normal(2)*(wplus+wminus)/2.0_rp

!    end function Riemann_Solver
!!----------------------------------------------------------------

!----------------------------------------------------------------
! Interpolate to Boundary
! Usage : real Interp_Value = Interpolate_to_Boundary(real phi, real l)
!----------------------------------------------------------------
    function Interpolate_to_Boundary(phi,l) result(Interp_Value)
        use Constants
        implicit none

        real(kind=rp), dimension(0:)            :: phi, l
        real(kind=rp)                           :: Interp_Value
        integer                                 :: N, j

        N = ubound(phi,1)
        Interp_Value = 0.0_rp
        do j = 0, N
            Interp_Value = Interp_Value + l(j)*phi(j)
        enddo
        
    end function Interpolate_to_Boundary
!----------------------------------------------------------------

!----------------------------------------------------------------
! A Matrix-Vector Multiplication Procedure
! Usage : MxVDerivative(real {D(i,j)}, complex {F}, int stind, int enind)
!----------------------------------------------------------------
! Variable Definition
! D      = Derivative Matrix
! F      = Function
!----------------------------------------------------------------

    function MxVDerivative(D, F) result(Fderiv)
        use Constants
        implicit none

        integer                                 :: N
        real(kind=rp), dimension(0:,0:)         :: D
        real(kind=rp), dimension(0:)            :: F
        real(kind=rp), dimension(0:ubound(F,1)) :: Fderiv
        integer                                 :: i, j
        real(kind=rp)                           :: t

        N = ubound(F,1)
        Fderiv(:) = 0.0_rp
        
        do i = 0, N
            t = 0.0_rp
            do j = 0, N
                t = t + D(i,j)*F(j)
            enddo
            Fderiv(i) = t
        enddo

    end function MxVDerivative
!----------------------------------------------------------------

end module NodalDG_2D_Class
