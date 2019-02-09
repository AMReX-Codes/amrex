module amrex_eb_flux_reg_nd_module
  use amrex_fort_module, only : amrex_real
  implicit none

  public

  integer, parameter :: crse_cell = 0
  integer, parameter :: crse_fine_boundary_cell = 1
  integer, parameter :: fine_cell = 2
  real(amrex_real), save :: reredistribution_threshold = 1.d-14

contains

  subroutine amrex_eb_disable_reredistribution () bind(c, name='amrex_eb_disable_reredistribution')
    reredistribution_threshold = 1.d10
  end subroutine amrex_eb_disable_reredistribution

end module amrex_eb_flux_reg_nd_module
