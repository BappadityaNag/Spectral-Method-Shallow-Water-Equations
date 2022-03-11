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
