module mg_tower_module

  use multifab_module
  use bl_constants_module

  implicit none

  integer, parameter :: MG_SMOOTHER_GS_RB  = 1
  integer, parameter :: MG_SMOOTHER_JACOBI = 2
  integer, parameter :: MG_SMOOTHER_MINION_CROSS = 5
  integer, parameter :: MG_SMOOTHER_MINION_FULL = 6
  integer, parameter :: MG_SMOOTHER_EFF_RB = 7

  integer, parameter :: MG_FCycle = 1
  integer, parameter :: MG_WCycle = 2
  integer, parameter :: MG_VCycle = 3
  integer, parameter :: MG_FVCycle = 4

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
     integer :: nub = 0
     integer :: cycle_type = MG_Vcycle

     ! bottom solver defaults good for bicg
     integer :: bottom_solver = 1
     integer :: bottom_max_iter = 100
     logical :: bottom_singular = .false.
     real(kind=dp_t) :: bottom_solver_eps = 1.0e-4_dp_t
     integer :: fancy_bottom_type = 1

     ! This must be true in order to enforce solvability
     logical :: coeffs_sum_to_zero = .false.

     integer :: ptype = 3
     logical :: use_lininterp = .true.

!     integer :: nboxes =  0
     integer :: nlevels =  0

     ! use hypre instead of mg_tower's multigrid for the solve itself
     integer :: use_hypre =  0

     ! allow code to subtract mean from RHS for singular problems
     logical :: ok_to_fix_singular = .true.

     ! let MG pick the maximum number of levels
     integer :: max_nlevel = 1024
     integer :: max_bottom_nlevel = 3
     integer :: min_width  = 2

     ! good for many problems
     integer :: max_iter = 50
     logical :: abort_on_max_iter = .true.
     logical :: always_use_bnorm  = .false.
     real(kind=dp_t) ::     eps = 1.0e-10_dp_t
     real(kind=dp_t) :: abs_eps = -ONE

     real(kind=dp_t) :: max_L0_growth = -ONE

     type(box), pointer :: pd(:) => Null()
     real(kind=dp_t), pointer :: dh(:,:) => Null()

     logical :: uniform_dh = .true.

     ! We now have pre-determined types of stencils with associated values of lcross
     integer :: stencil_type

     ! If lcross is true and passed into multifab_fill_boundary, then 
     !    multifab_fill_boundary knows not to fill any corner cells.
     logical :: lcross

     type(multifab),  pointer :: cc(:) => Null()
     type(multifab),  pointer :: ff(:) => Null()
     type(multifab),  pointer :: dd(:) => Null()
     type(multifab),  pointer :: uu(:) => Null()
     type(multifab),  pointer :: ss(:) => Null()
     type(imultifab), pointer :: mm(:) => Null()

     integer, pointer :: face_type(:,:,:)  => Null()
     logical, pointer :: skewed(:,:)       => Null()
     logical, pointer :: skewed_not_set(:) => Null()
     integer, pointer :: domain_bc(:,:)    => Null()

     ! Only relevant if bottom_solver == 1, 2 or 3 AND nodal
     type(multifab) :: nodal_mask

     integer ::    verbose = 0
     integer :: cg_verbose = 0

     type(mg_tower), pointer :: bottom_mgt  => Null()
     integer,        pointer :: bottom_comm => Null()

  end type mg_tower

end module mg_tower_module
