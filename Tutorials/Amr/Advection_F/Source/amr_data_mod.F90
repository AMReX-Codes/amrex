
module amr_data_module

  use iso_c_binding
  use amrex_amr_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_particlecontainer_module, only: amrex_particlecontainer, &
       amrex_particlecontainer_destroy
  
  implicit none

  private
  public :: t_new, t_old, phi_new, phi_old, flux_reg, pc
  public :: amr_data_init, amr_data_finalize

  real(rt), allocatable :: t_new(:)
  real(rt), allocatable :: t_old(:)

  type(amrex_multifab), allocatable :: phi_new(:)
  type(amrex_multifab), allocatable :: phi_old(:)

  type(amrex_fluxregister), allocatable :: flux_reg(:)

  type(amrex_particlecontainer) :: pc
  
contains

  subroutine amr_data_init ()
    allocate(t_new(0:amrex_max_level))
    t_new = 0.0_rt

    allocate(t_old(0:amrex_max_level))
    t_old = -1.0e100_rt

    allocate(phi_new(0:amrex_max_level))
    allocate(phi_old(0:amrex_max_level))

    allocate(flux_reg(0:amrex_max_level))
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
    call amrex_particlecontainer_destroy(pc)
  end subroutine amr_data_finalize
  
end module amr_data_module
