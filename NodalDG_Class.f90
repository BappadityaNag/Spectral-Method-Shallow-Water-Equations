!----------------------------------------------------------------
! Contains 
! 1. Construct
! 2. Destruct
! 3. DG_Space_Derivative
! 4. DG_Time_Derivative
! 5. Interpolate_to_Boundary
! 6. DGStep_by_RK3
! 7. Nodal_DG_Integrator
! 8. MxVDerivative
!----------------------------------------------------------------

module NodalDG_Class
implicit none
contains

!----------------------------------------------------------------
! Construct the Contour
! Usage : gamma = Construct(nEqn, N, M)
!----------------------------------------------------------------
    function Construct(N) result(gamma)
        use Constants
        use Quadrature_Points_and_Weights_Class
        use Lagrange_Interpolation_Class
        implicit none

        integer                                    :: N
        type(NodalDG_Storage)                      :: gamma
        integer                                    :: i, j
        real(kind=rp), dimension(:,:), allocatable :: D
        real(kind=rp), dimension(0:N)              :: WB
        
        gamma%N = N
        allocate(gamma%phi(0:N))
        allocate(gamma%X(0:N))
        allocate(gamma%W(0:N))
        allocate(gamma%D1(0:N,0:N))
        allocate(gamma%D2(0:N,0:N))
        allocate(gamma%lleft(0:N))
        allocate(gamma%lright(0:N))
        allocate(D(0:N,0:N))
        
        call Legendre_Gauss_Nodes_and_Weights(N,gamma%X,gamma%W)
        WB = Barycentric_Weights(gamma%X)
        gamma%lleft = Lagrange_Interpolating_Polynomials(-1.0_rp,gamma%X,WB)
        gamma%lright = Lagrange_Interpolating_Polynomials(1.0_rp,gamma%X,WB)
        D = Polynomial_Derivative_Matrix(gamma%X)
        
        do j = 0,N
            do i = 0,N
                gamma%D1(i,j) = -D(j,i)*gamma%W(j)/gamma%W(i)
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

        type(NodalDG_Storage)                      :: gamma

        deallocate(gamma%phi)
        deallocate(gamma%X)
        deallocate(gamma%W)
        deallocate(gamma%D1)
        deallocate(gamma%D2)
        deallocate(gamma%lleft)
        deallocate(gamma%lright)

    end subroutine Destruct
!----------------------------------------------------------------

!----------------------------------------------------------------
! First Spatial Derivative
! Usage : {phiDeriv} = DG_Space_Derivative({phiL},{phiR},gamma)
!----------------------------------------------------------------
    function DG_Space_Derivative(phiL,phiR,gamma) result(phiDeriv)
        use Constants
        implicit none

        real(kind=rp)                                     :: phiL, phiR
        type(NodalDG_Storage)                             :: gamma
        real(kind=rp), dimension(0:ubound(gamma%phi,1))   :: phiDeriv
        integer                                           :: N, j

        N     = ubound(gamma%phi,1)
        phiDeriv = MxVDerivative(gamma%D1,gamma%phi)
        do j = 0, N
            phiDeriv(j) = phiDeriv(j) + ( phiR*gamma%lright(j) - phiL*gamma%lleft(j) ) / gamma%W(j)
        enddo

    end function DG_Space_Derivative
!----------------------------------------------------------------

!----------------------------------------------------------------
! Time Derivative
! Usage : {Qdot} = DG_Time_Derivative(t)
!----------------------------------------------------------------
    function DG_Time_Derivative(t,gamma,c,bc) result(phiDot)
        use Constants
        implicit none

        real(kind=rp)                                     :: t, c, bc, phiL, phiR
        type(NodalDG_Storage)                             :: gamma
        real(kind=rp), dimension(0:ubound(gamma%phi,1))   :: phiDot

        if (c .gt. 0.0_rp) then
            phiL = bc
            phiR = Interpolate_to_Boundary(gamma%phi,gamma%lright)
        else
            phiR = bc
            phiL = Interpolate_to_Boundary(gamma%phi,gamma%lleft)
        endif
        phiDot = -c * DG_Space_Derivative(phiL,phiR,gamma)

    end function DG_Time_Derivative
!----------------------------------------------------------------

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
! Low Storage Runge-Kutta Integration of a Nodal DG Approximation
! Usage : DGStep_by_RK3(P,D,tn,dt)
!----------------------------------------------------------------
    subroutine DGStep_by_RK3(tn,dt,gamma,c,bc)
        use Constants
        use Lagrange_Interpolation_Class
        implicit none

        type(NodalDG_Storage)                           :: gamma
        real(kind=rp), dimension(0:ubound(gamma%phi,1)) :: phiDot, tmp
        real(kind=rp)                                   :: t, tn, dt
        integer                                         :: N, j, k
        real(kind=rp)                                   :: c, bc

        real(kind=rp), dimension(3) :: a = (/ 0.0_rp, -5.0_rp/9.0_rp, -153.0_rp/128.0_rp /)
        real(kind=rp), dimension(3) :: b = (/ 0.0_rp, 1.0_rp/3.0_rp, 3.0_rp/4.0_rp /)
        real(kind=rp), dimension(3) :: g = (/ 1.0_rp/3.0_rp, 15.0_rp/16.0_rp, 8.0_rp/15.0_rp /)

        N    = ubound(gamma%phi,1)
        tmp = 0.0_rp
        do k = 1,3
            t = tn + b(k)*dt
            phiDot = DG_Time_Derivative(t,gamma,c,bc)
            do j = 0,N
                tmp(j) = a(k)*tmp(j) + phiDot(j)
                gamma%phi(j) = gamma%phi(j) + g(k)*dt*tmp(j)
            enddo
        enddo

    end subroutine DGStep_by_RK3
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

        integer                                 :: N, i, j
        real(kind=rp), dimension(0:,0:)         :: D
        real(kind=rp), dimension(0:)            :: F
        real(kind=rp), dimension(0:ubound(F,1)) :: Fderiv
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

end module NodalDG_Class
