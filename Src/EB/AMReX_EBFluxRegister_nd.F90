#include <AMReX_Config.H>

module amrex_eb_flux_reg_nd_module
  use amrex_fort_module, only : amrex_real
  implicit none

  public

  integer, parameter :: crse_cell = 0
  integer, parameter :: crse_fine_boundary_cell = 1
  integer, parameter :: fine_cell = 2
#ifdef AMREX_USE_FLOAT
  real(amrex_real), save :: reredistribution_threshold = 1.e-5_amrex_real
#else
  real(amrex_real), save :: reredistribution_threshold = 1.d-14
#endif

contains

  subroutine amrex_eb_disable_reredistribution () bind(c, name='amrex_eb_disable_reredistribution')
    reredistribution_threshold = 1.d10
  end subroutine amrex_eb_disable_reredistribution

  real(amrex_real) function amrex_eb_get_reredistribution_threshold () &
       bind(c, name='amrex_eb_get_reredistribution_threshold')
    amrex_eb_get_reredistribution_threshold = reredistribution_threshold
  end function amrex_eb_get_reredistribution_threshold

end module amrex_eb_flux_reg_nd_module
