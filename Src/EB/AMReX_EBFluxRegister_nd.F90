module amrex_eb_flux_reg_nd_module
  use amrex_fort_module, only : amrex_real
  use amrex_ya_flux_reg_nd_module, only : crse_cell, crse_fine_boundary_cell, fine_cell
  implicit none

  public

  real(amrex_real), save, public :: reredistribution_threshold = 1.d-14

contains

  subroutine amrex_eb_disable_reredistribution () bind(c, name='amrex_eb_disable_reredistribution')
    reredistribution_threshold = 1.d10
  end subroutine amrex_eb_disable_reredistribution

end module amrex_eb_flux_reg_nd_module
