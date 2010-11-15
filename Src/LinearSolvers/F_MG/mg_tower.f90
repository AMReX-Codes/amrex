module mg_tower_module

  use multifab_module
  use stencil_module
! use sparse_solve_module
  use bl_timer_module

  implicit none

  integer, parameter :: MG_SMOOTHER_GS_RB  = 1
  integer, parameter :: MG_SMOOTHER_JACOBI = 2
  integer, parameter :: MG_SMOOTHER_GS_LEX = 3
  integer, parameter :: MG_SMOOTHER_MINION_CROSS = 5
  integer, parameter :: MG_SMOOTHER_MINION_FULL = 6
  integer, parameter :: MG_SMOOTHER_EFF_RB = 7

  integer, parameter :: MG_FCycle = 1
  integer, parameter :: MG_WCycle = 2
  integer, parameter :: MG_VCycle = 3

  type mg_tower

     integer :: dim = 0

     ! defaults appropriate for abec, laplacian
     integer :: nc = 1
     integer :: ng = 1
     integer :: ns = 1

     ! gsrb
     integer :: smoother = MG_SMOOTHER_GS_RB
     integer :: nu1 = 2
     integer :: nu2 = 2
     integer :: nuf = 8
     integer :: nub = 10
     integer :: gamma = 1
     integer :: cycle_type = MG_Vcycle
     real(kind=dp_t) :: omega = 1.0_dp_t

     ! bottom solver defaults good for bicg
     integer :: bottom_solver = 1
     integer :: bottom_max_iter = 100
     logical :: bottom_singular = .false.
     real(kind=dp_t) :: bottom_solver_eps = 1.0e-4_dp_t

     integer :: nboxes =  0
     integer :: nlevels =  0

     ! use hypre instead of mg_tower's multigrid for the solve itself
     integer :: use_hypre =  0

     ! let MG pick the maximum number of levels
     integer :: max_nlevel = 1024
     integer :: max_bottom_nlevel = 3
     integer :: min_width  = 2

     ! good for many problems
     integer :: max_iter = 50
     real(kind=dp_t) ::     eps = 1.0e-10_dp_t
     real(kind=dp_t) :: abs_eps = -1.0_dp_t

     type(box), pointer :: pd(:) => Null()
     real(kind=dp_t), pointer :: dh(:,:) => Null()
     logical :: uniform_dh = .true.

     type(multifab), pointer :: cc(:) => Null()
     type(multifab), pointer :: ff(:) => Null()
     type(multifab), pointer :: dd(:) => Null()
     type(multifab), pointer :: uu(:) => Null()
     type(multifab), pointer :: ss(:) => Null()
     type(imultifab), pointer :: mm(:) => Null()
     type(stencil), pointer :: st(:) => Null()

     integer, pointer :: face_type(:,:,:)  => Null()
     logical, pointer :: skewed(:,:)       => Null()
     logical, pointer :: skewed_not_set(:) => Null()

     type(timer), pointer :: tm(:) => Null()

     ! Only relevant if bottom_solver == 3
!    type(sparse) sparse_object
!    type(multifab) :: rh1
!    type(multifab) :: uu1
!    type(multifab) :: ss1
!    type(imultifab) :: mm1

     ! Only relevant if bottom_solver == 1 or 2 AND nodal
     type(multifab) :: nodal_mask

     integer ::    verbose = 0
     integer :: cg_verbose = 0
     integer ::    st_type = 0

      type(mg_tower), pointer :: bottom_mgt => Null()

  end type mg_tower

end module mg_tower_module
