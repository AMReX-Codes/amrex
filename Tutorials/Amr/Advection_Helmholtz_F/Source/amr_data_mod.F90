
module amr_data_module

  use iso_c_binding
  use amrex_amr_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_base_module
  use amrex_linear_solver_module
 
  implicit none

  logical, save :: composite_solve = .true.

  integer, save :: verbose = 2
  integer, save :: cg_verbose = 0
  integer, save :: max_iter = 100
  integer, save :: max_fmg_iter = 0
  integer, save :: linop_maxorder = 2
  logical, save :: agglomeration = .true.
  logical, save :: consolidation = .true.
  real(rt), allocatable :: t_new(:)
  real(rt), allocatable :: t_old(:)

  type(amrex_multifab), allocatable :: phi_new(:)
  type(amrex_multifab), allocatable :: phi_old(:)
  
  type(amrex_fluxregister), allocatable :: flux_reg(:)

  type(amrex_boxarray), allocatable, save :: ba_new(:)
  type(amrex_distromap), allocatable, save :: dm_new(:)

  type(amrex_multifab), allocatable, save :: solution(:)
  type(amrex_multifab), allocatable, save :: rhs(:)
  type(amrex_multifab), allocatable, save :: exact_solution(:)
  type(amrex_multifab), allocatable, save :: acoef(:)
  type(amrex_multifab), allocatable, save :: bcoef(:)
  type(amrex_multifab),allocatable :: beta(:,:)
  real(amrex_real), save :: ascalar, bscalar
  type(amrex_abeclaplacian) :: abeclap

  integer :: ncomp

contains

  subroutine amr_data_init ()

    allocate(t_new(0:amrex_max_level))
    t_new = 0.0_rt

    allocate(t_old(0:amrex_max_level))
    t_old = -1.0e100_rt

    allocate(phi_new(0:amrex_max_level))
    allocate(phi_old(0:amrex_max_level))

    allocate(flux_reg(0:amrex_max_level))
    
    allocate(ba_new(0:amrex_max_level))
    allocate(dm_new(0:amrex_max_level))

    allocate(solution(0:amrex_max_level))
    allocate(rhs(0:amrex_max_level))
    allocate(exact_solution(0:amrex_max_level))
    allocate(acoef(0:amrex_max_level))
    allocate(bcoef(0:amrex_max_level))
    allocate(beta(amrex_spacedim,0:amrex_max_level))

  end subroutine amr_data_init

  subroutine amr_data_finalize
    integer :: lev
    do lev = 0, amrex_max_level
       call amrex_multifab_destroy(phi_new(lev))
       call amrex_multifab_destroy(phi_old(lev))
    end do
    do lev = 1, amrex_max_level
       call amrex_fluxregister_destroy(flux_reg(lev))
    end do
  end subroutine amr_data_finalize
  
end module amr_data_module
