!----------------------------------------------------------------
! Contains
!  1. Almost_Equal
!  2. Barycentric_Weights
!  3. Lagrange_Interpolation
!  4. Polynomial_Interpolation_Matrix
!  5. Interpolate_to_New_Points
!  6. Lagrange_Interpolating_Polynomials
!  7. Lagrange_Interpolant_Derivative
!  8. Polynomial_Derivative_Matrix
!  9. mth_Order_Polynomial_Derivative_Matrix
! 10. MxV_Derivative
!----------------------------------------------------------------

module Lagrange_Interpolation_Class
implicit none
contains

!----------------------------------------------------------------
! Testing Equality of Two Floating Point Numbers
! Usage : true/false = Almost_Equal(x,y)
!----------------------------------------------------------------

    function Almost_Equal(x, y)
        use Constants
        implicit none

        real(kind=rp)                           :: x, y
        logical                                 :: Almost_Equal

        if ( x .eq. 0.0_rp .or. y .eq. 0.0_rp ) then
            if ( abs(x-y) .le. 2.0_rp*epsilon(x) ) then
                Almost_Equal = .True.
            else
                Almost_Equal = .False.
            endif
        else
            if ( abs(x-y) .le. epsilon(x)*abs(x) .and. abs(x-y) .le. epsilon(x)*abs(y) ) then
                Almost_Equal = .True.
            else
                Almost_Equal = .False.
            endif
        endif
        
    end function Almost_Equal
!----------------------------------------------------------------

!----------------------------------------------------------------
! Weights for Lagrange Interpolation
! Usage : real {W} = Barycentric_Weights(real {X})
!----------------------------------------------------------------

    function Barycentric_Weights(X) result(W)
        use Constants
        implicit none

        real(kind=rp), dimension(0:)            :: X
        real(kind=rp), dimension(0:ubound(X,1)) :: W
        integer                                 :: N, j, k

        N = ubound(X,1)
        W(:) = 1.0_rp
        do j = 1, N
            do k = 0, j-1
                W(k) = W(k)*(X(k)-X(j))
                W(j) = W(j)*(X(j)-X(k))
            enddo
        enddo
        do j = 0, N
            W(j) = 1.0_rp/W(j)
        enddo

    end function Barycentric_Weights
!----------------------------------------------------------------

!----------------------------------------------------------------
! Lagrange Interpolant from Barycentric Form
! Usage : interp = Lagrange_Interpolation(real x, real {X}, real {f}, real {W})
!----------------------------------------------------------------

    function Lagrange_Interpolation(x, Xarr, f, W) result(interp)
        use Constants
        implicit none

        real(kind=rp), dimension(0:)            :: Xarr, f, W
        integer                                 :: N, j
        real(kind=rp)                           :: x, t, numer, denom, interp

        N = ubound(Xarr,1)
        numer = 0.0_rp
        denom = 0.0_rp
        do j = 0, N
            if (Almost_Equal(x,Xarr(j))) then
                interp = f(j)
                return
            endif
            t = W(j)/(x-Xarr(j))        
            numer = numer + t*f(j)
            denom = denom + t
        enddo
        interp = numer/denom
        
    end function Lagrange_Interpolation
!----------------------------------------------------------------

!----------------------------------------------------------------
! Matrix for Interpolation Between Two Sets of Points
! Usage : real {T} = Polynomial_Interpolation_Matrix(real {X}, real {W}, real {xi})
!----------------------------------------------------------------

    function Polynomial_Interpolation_Matrix(X,W,xi) result(T)
        use Constants
        implicit none

        real(kind=rp), dimension(0:)            :: X, W, xi
        real(kind=rp), dimension(0:ubound(xi,1),0:ubound(X,1)) :: T
        integer                                 :: M, N, j, k
        real(kind=rp)                           :: s, u
        logical                                 :: rowHasMatch

        M = ubound(xi,1)
        N = ubound(X,1)
        do k = 0, M
            rowHasMatch = .False.
            do j = 0, N
                T(k,j) = 0.0_rp
                if (Almost_Equal(xi(k),X(j))) then
                    rowHasMatch = .True.
                    T(k,j) = 1.0_rp
                endif
            enddo
            if (rowHasMatch .eqv. .False.) then
                s = 0.0_rp
                do j = 0, N
                    u = W(j)/(xi(k)-X(j))
                    T(k,j) = u
                    s = s + u
                enddo
                do j = 0, N
                    T(k,j) = T(k,j)/s
                enddo
            endif
        enddo
        
    end function Polynomial_Interpolation_Matrix
!----------------------------------------------------------------

!----------------------------------------------------------------
! Interpolation between Two sets of Points by Matrix Multiplication
! Usage : {fInterp} = Interpolate_to_New_Points(real {T}, real {f})
!----------------------------------------------------------------

    function Interpolate_to_New_Points(T,f) result(fInterp)
        use Constants
        implicit none

        real(kind=rp), dimension(0:,0:)         :: T
        real(kind=rp), dimension(0:)            :: f
        real(kind=rp), dimension(0:ubound(T,1)) :: fInterp
        integer                                 :: M, N, j, k
        real(kind=rp)                           :: u

        M = ubound(T,1)
        N = ubound(f,1)
        do k = 0, M
            u = 0.0_rp
            do j = 0, N
                u = u + T(k,j) * f(j)
            enddo
            fInterp(k) = u
        enddo

    end function Interpolate_to_New_Points
!----------------------------------------------------------------

!----------------------------------------------------------------
! Lagrange Interpolating Polynomials : lj(x)
! Usage : {l} = Lagrange_Interpolating_Polynomials(real x, real {Xarr}, real {W})
!----------------------------------------------------------------

    function Lagrange_Interpolating_Polynomials(x,Xarr,W) result(l)
        use Constants
        implicit none

        real(kind=rp), dimension(0:)               :: Xarr, W
        real(kind=rp), dimension(0:ubound(Xarr,1)) :: l
        real(kind=rp)                              :: x, s, t
        integer                                    :: N, j
        logical                                    :: xMatchesNode

        N = ubound(Xarr,1)
        xMatchesNode = .False.
        do j = 0, N
            l(j) = 0.0_rp
            if ( Almost_Equal( x, Xarr(j) )) then
                l(j) = 1.0_rp
                xMatchesNode = .True.
            endif
        enddo
        if (xMatchesNode) return
        s = 0.0_rp
        do j = 0, N
            t = W(j) / ( x - Xarr(j) )
            l(j) = t
            s = s + t
        enddo
        do j = 0, N
            l(j) = l(j)/s
        enddo
            
    end function Lagrange_Interpolating_Polynomials
!----------------------------------------------------------------

!----------------------------------------------------------------
! Direct Computation of the Polynomial Derivative in Barycentric Form
! Usage : {pderiv} = Lagrange_Interpolant_Derivative(real x, real {Xarr}, real {f}, real {W})
!----------------------------------------------------------------

    function Lagrange_Interpolant_Derivative(x,Xarr,f,W) result(pderiv)
        use Constants
        implicit none

        real(kind=rp), dimension(0:)               :: Xarr, f, W
        real(kind=rp)                              :: x, p, pderiv, numer, denom, t
        integer                                    :: N, i, j
        logical                                    :: atNode

        N = ubound(Xarr,1)
        atNode = .False.
        numer = 0.0_rp
        do j = 0, N
            if ( Almost_Equal( x, Xarr(j) )) then
                atNode = .True.
                p = f(j)
                denom = -W(j)
                i = j
            endif
        enddo
        if (atNode) then
            do j = 0, N
                if (j .ne. i) numer = numer + W(j) * ( p - f(j) ) / ( x - Xarr(j) )
            enddo
        else
            denom = 0.0_rp
            p = Lagrange_Interpolation(x, Xarr, f, W)
            do j = 0, N
                t = W(j) / ( x - Xarr(j) )
                numer = numer + t * ( p - f(j) ) / ( x - Xarr(j) )
                denom = denom + t
            enddo
        endif
        pderiv = numer/denom

    end function Lagrange_Interpolant_Derivative
!----------------------------------------------------------------

!----------------------------------------------------------------
! First Order Polynomial Derivative Matrix
! Usage : {D} = Polynomial_Derivative_Matrix(real {X})
!----------------------------------------------------------------

    function Polynomial_Derivative_Matrix(X) result(D)
        use Constants
        implicit none

        real(kind=rp), dimension(0:)               :: X
        real(kind=rp), dimension(0:ubound(X,1),0:ubound(X,1)) :: D
        real(kind=rp), dimension(0:ubound(X,1))    :: W
        integer                                    :: N, i, j

        N = ubound(X,1)
        W = Barycentric_Weights(X)
        do i = 0, N
            D(i,i) = 0.0_rp
            do j = 0, N
                if (j .ne. i) then
                    D(i,j) = W(j) / ( W(i) * ( X(i) - X(j) ) )
                    D(i,i) = D(i,i) - D(i,j)
                endif
            enddo
        enddo

    end function Polynomial_Derivative_Matrix
!----------------------------------------------------------------

!----------------------------------------------------------------
! mth Order Polynomial Derivative Matrix
! Usage : {D} = mth_Order_Polynomial_Derivative_Matrix(int m, real {X})
!----------------------------------------------------------------

    function mth_Order_Polynomial_Derivative_Matrix(m,X) result(D)
        use Constants
        implicit none

        integer                                    :: N, m, i, j, k
        real(kind=rp), dimension(0:)               :: X
        real(kind=rp), dimension(0:ubound(X,1),0:ubound(X,1)) :: D
        real(kind=rp), dimension(0:ubound(X,1))    :: W, tmp

        N = ubound(X,1)
        W = Barycentric_Weights(X)
        D = Polynomial_Derivative_Matrix(X)
        
        if (m .eq. 1) return
        do k = 2, m
            do i = 0,N
                tmp(i) = D(i,i)
                D(i,i) = 0.0_rp
                do j = 0, N
                    if (j .ne. i) then
                        D(i,j) = real(k,kind=rp)/(X(i)-X(j)) * (W(j)/W(i)*tmp(i)-D(i,j))
                        D(i,i) = D(i,i) - D(i,j)
                    endif
                enddo
            enddo
        enddo

    end function mth_Order_Polynomial_Derivative_Matrix
!----------------------------------------------------------------

!----------------------------------------------------------------
! A Matrix-Vector Multiplication Procedure
! Usage : real {F'} = MxV_Derivative(real {D}, real {F})
!----------------------------------------------------------------

    function MxV_Derivative(D, F) result(Fderiv)
        use Constants
        implicit none

        real(kind=rp), dimension(0:,0:)         :: D
        real(kind=rp), dimension(0:)            :: F
        real(kind=rp), dimension(0:ubound(F,1)) :: Fderiv
        integer                                 :: N, i, j
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

    end function MxV_Derivative
!----------------------------------------------------------------

end module Lagrange_Interpolation_Class
