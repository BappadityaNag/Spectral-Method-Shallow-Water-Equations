!----------------------------------------------------------------
! Defines Precision and Constants
!----------------------------------------------------------------

module Constants

    implicit none
    
    ! Precision
    integer, parameter           :: rp = selected_real_kind(15)

    ! Parameters
    integer, parameter           :: MaxIter = 100000000
    real(kind=rp), parameter     :: pi = acos(-1.0_rp)
    real(kind=rp),parameter      :: tol = 1.0e-10
    integer, parameter           :: Null = -999

    type :: grid
        real                                              :: delta
        real(kind=rp), dimension(:,:), allocatable        :: G
    end type grid

    type :: edge
        integer                                           :: N
        real(kind=rp), dimension(:), allocatable          :: nodes, X, Y, W
    end type edge

    type :: NodalDG_Storage
        integer                                           :: N
        real(kind=rp), dimension(:), allocatable          :: X, W, phi
        real(kind=rp), dimension(:,:), allocatable        :: D1, D2
        real(kind=rp), dimension(:), allocatable          :: lleft, lright
    end type NodalDG_Storage

    type :: NodalDG_2DStorage
        integer                                           :: N, M
        real(kind=rp), dimension(:), allocatable          :: xi, Wxi, WBxi
        real(kind=rp), dimension(:,:), allocatable        :: D1xi, D2xi
        real(kind=rp), dimension(:), allocatable          :: eta, Weta, WBeta
        real(kind=rp), dimension(:,:), allocatable        :: D1eta, D2eta
        real(kind=rp), dimension(:), allocatable          :: lxileft, lxiright, letaleft, letaright
    end type NodalDG_2DStorage

    type :: NodalDG_2DSystem
        integer                                           :: nEqn
        type(NodalDG_2DStorage)                           :: spA
        real(kind=rp), dimension(:,:,:), allocatable      :: Q
    end type NodalDG_2DSystem

! Shared_Node_Ptrs used in SEM1D_Class
    type :: Shared_Node_Ptrs
        integer                                           :: eLeft, eRight, nodeLeft, nodeRight
    end type Shared_Node_Ptrs

! SEM1D Element
! Data Storage for the 1D-SEM
    type :: SEM1D
        integer                                           :: N, K
        real(kind=rp), dimension(:), allocatable          :: x              ! Element Boundary Locations
        real(kind=rp), dimension(:), allocatable          :: xi             ! Node Locations
        real(kind=rp), dimension(:), allocatable          :: W              ! Legendre Gauss-Lobatto quadrature weights
        real(kind=rp), dimension(:), allocatable          :: deltax         ! Element Sizes
        real(kind=rp), dimension(:,:), allocatable        :: G              ! Derivative Matrix
        type(Shared_Node_Ptrs), dimension(:), allocatable :: p              ! Shared Node Pointers
        real(kind=rp), dimension(:,:), allocatable        :: phi            ! Solution
        real(kind=rp), dimension(:,:), allocatable        :: RHS            ! For implicit Integration
    end type SEM1D

! Element in 1D
    type :: Element
        integer                                           :: N
        integer                                           :: nEqn           ! # of eqns. in the physical system
        real(kind=rp)                                     :: deltax, xL, xR ! Size and left and right boundaries of the element
        type(NodalDG_Storage)                             :: dG
        real(kind=rp), dimension(:,:), allocatable        :: phi, phidot    ! Solution and time derivative
        real(kind=rp), dimension(:,:), allocatable        :: G              ! For low-storage Runga-Kutta
        real(kind=rp), dimension(:), allocatable          :: phiL, phiR     ! Solution on left/right element boundary
        real(kind=rp), dimension(:), allocatable          :: FL, FR         ! Flux on left/right element boundary
    end type Element

! Mesh 1D
    type :: Mesh1D
        integer                                           :: K              ! # of Elements
        type(Element), dimension(:), allocatable          :: e              ! Elements
        type(Shared_Node_Ptrs), dimension(:), allocatable :: p              ! Shared Node Pointers
    end type Mesh1D

end module Constants
