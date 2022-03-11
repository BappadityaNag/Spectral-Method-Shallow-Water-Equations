!----------------------------------------------------------------
! Contains 
! 1. qAndL_Evaluation
! 2. Legendre_GaussLobatto_Nodes_and_Weights
! 3. Legendre_Polynomial_and_Derivative
! 4. Legendre_Gauss_Nodes_and_Weights
!----------------------------------------------------------------

module Quadrature_Points_and_Weights_Class
implicit none
contains

!----------------------------------------------------------------
! Combined Algorithm to Compute LN(x), q(x) = L(N+1) − L(N−1), and q'(x)
! Usage : qAndL_Evaluation(int N, real x, real q, real qderiv, real LN)
!----------------------------------------------------------------

    subroutine qAndL_Evaluation(N,x,q,qderiv,LN)
        use Constants
        implicit none

        integer                                :: N, k
        real(kind=rp)                          :: x, q, qderiv, LN
        real(kind=rp), dimension(N-2:N+1)      :: L, Lderiv

        L(N-2) = 1.0_rp
        L(N-1) = x
        Lderiv(N-2) = 0.0_rp
        Lderiv(N-1) = 1.0_rp
        
        do k = 2, N
            L(N) = (2.0_rp*k-1.0_rp)/k * x*L(N-1) - (real(k,kind=rp)-1.0_rp)/k * L(N-2)
            Lderiv(N) = Lderiv(N-2) + (2.0_rp*k-1.0_rp)*L(N-1)
            L(N-2) = L(N-1)
            L(N-1) = L(N)
            Lderiv(N-2) = Lderiv(N-1)
            Lderiv(N-1) = Lderiv(N)
        enddo
        
        k = N + 1
        L(N+1) = (2.0_rp*k-1.0_rp)/k * x*L(N) - (1.0_rp*k-1.0_rp)/k * L(N-2)
        Lderiv(N+1) = Lderiv(N-2) + (2.0_rp*k-1.0_rp)*L(N-1)
        q = L(N+1) - L(N-2)
        qderiv = Lderiv(N+1) - Lderiv(N-2)
        LN = L(N)
        
    end subroutine qAndL_Evaluation
!----------------------------------------------------------------

!----------------------------------------------------------------
! Usage : Legendre_GaussLobatto_Nodes_and_Weights(int N, real {X}, real {W})
!----------------------------------------------------------------

    subroutine Legendre_GaussLobatto_Nodes_and_Weights(N,X,W)
        use Constants
        implicit none

        integer                                :: N, j, k
        real(kind=rp), dimension(0:N)          :: X, W
        real(kind=rp)                          :: q, qderiv, LN, del
        
        if (N .eq. 1) then
            X(0) = -1.0_rp
            W(0) = 1.0_rp
            X(1) = 1.0_rp
            W(1) = W(0)
        else
            X(0) = -1.0_rp
            W(0) = 2.0_rp / (N*(N+1))
            X(N) = 1.0_rp
            W(N) = W(0)
            do j = 1, floor( (real(N,kind=rp)+1.0_rp)/2.0_rp) - 1
                X(j) = -cos( (j+0.25_rp)*pi/N - 3.0_rp/(8.0_rp*N*pi*(j+0.25_rp)) )
                do k = 0, MaxIter
                    call qAndL_Evaluation(N,X(j),q,qderiv,LN)
                    del = -q/qderiv
                    X(j) = X(j) + del
                    if (abs(del) .le. epsilon(X(j))*abs(X(j))) exit
                enddo
                call qAndL_Evaluation(N,X(j),q,qderiv,LN)
                X(N-j) = -X(j)
                W(j) = 2.0_rp/(N*(N+1)*LN**2)
                W(N-j) = W(j)
            enddo
        endif
        if (mod(N,2) .eq. 0) then
            call qAndL_Evaluation(N,0.0_rp,q,qderiv,LN)
            X(N/2) = 0.0_rp
            W(N/2) = 2.0_rp/(N*(N+1)*LN**2)
        endif

    end subroutine Legendre_GaussLobatto_Nodes_and_Weights
!----------------------------------------------------------------

!----------------------------------------------------------------
! The Legendre Polynomial of Degree N and its Derivative using the Three
! Term  Recursion
! Usage : Legendre_Polynomial_and_Derivative(int N, real x, real LN, real LNderiv)
!----------------------------------------------------------------

    subroutine Legendre_Polynomial_and_Derivative(N,x,LN,LNderiv)
        use Constants
        implicit none

        integer                                :: N, k
        real(kind=rp)                          :: x, LN, LNderiv
        real(kind=rp), dimension(N-2:N)        :: L, Lderiv

        if (N .eq. 0) then
            L(N)        = 1.0_rp
            Lderiv(N)   = 0.0_rp
        elseif (N .eq. 1) then
            L(N)        = x
            Lderiv(N)   = 1.0_rp
        else
            L(N-2)      = 1.0_rp
            L(N-1)      = x
            Lderiv(N-2) = 0.0_rp
            Lderiv(N-1) = 1.0_rp
            do k = 2, N
                L(N) = (2.0_rp*k-1.0_rp)/k * x*L(N-1) - (real(k,kind=rp)-1.0_rp)/k * L(N-2)
                Lderiv(N) = Lderiv(N-2) + (2.0_rp*k-1.0_rp)*L(N-1)
                L(N-2) = L(N-1)
                L(N-1) = L(N)
                Lderiv(N-2) = Lderiv(N-1)
                Lderiv(N-1) = Lderiv(N)
            enddo
        endif        
        LN = L(N)
        LNderiv = Lderiv(N)
        
    end subroutine Legendre_Polynomial_and_Derivative
!----------------------------------------------------------------

!----------------------------------------------------------------
! Usage : Legendre_Gauss_Nodes_and_Weights(int N, real {X}, real {W})
!----------------------------------------------------------------

    subroutine Legendre_Gauss_Nodes_and_Weights(N,X,W)
        use Constants
        implicit none

        integer                         :: N, j, k
        real(kind=rp), dimension(0:N)   :: X, W
        real(kind=rp)                   :: LNplus1, LNplus1deriv, del
        
        if (N .eq. 0) then
            X(0) = 0.0_rp
            W(0) = 2.0_rp
        elseif (N .eq. 1) then
            X(0) = -sqrt(1.0_rp/3.0_rp)
            W(0) = 1.0_rp
            X(1) = -X(0)
            W(1) = W(0)
        else
            do j = 0, floor( (real(N,kind=rp)+1.0_rp)/2.0_rp ) - 1
                X(j) = -cos( (2.0_rp*j+1.0_rp)*pi/(2.0_rp*N+2.0_rp) )
                do k = 0, MaxIter
                    call Legendre_Polynomial_and_Derivative(N+1,X(j),LNplus1,LNplus1deriv)
                    del = -LNplus1/LNplus1deriv
                    X(j) = X(j) + del
                    if (abs(del) .le. epsilon(X(j))*abs(X(j))) exit
                enddo
                call Legendre_Polynomial_and_Derivative(N+1,X(j),LNplus1,LNplus1deriv)
                X(N-j) = -X(j)
                W(j) = 2.0_rp/( (1.0_rp-X(j)**2)*(LNplus1deriv)**2 )
                W(N-j) = W(j)
            enddo
        endif
        if (mod(N,2) .eq. 0) then
            call Legendre_Polynomial_and_Derivative(N+1,0.0_rp,LNplus1,LNplus1deriv)
            X(N/2) = 0.0_rp
            W(N/2) = 2.0_rp/(LNplus1deriv**2)
        endif

    end subroutine Legendre_Gauss_Nodes_and_Weights
!----------------------------------------------------------------

end module Quadrature_Points_and_Weights_Class
