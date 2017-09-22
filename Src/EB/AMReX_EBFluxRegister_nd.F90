module amrex_eb_flux_reg_nd_module
  use amrex_fort_module, only : amrex_real
  use amrex_ya_flux_reg_nd_module, only : crse_cell, crse_fine_boundary_cell, fine_cell
  implicit none

  public

  real(amrex_real), parameter, public :: reredistribution_threshold = 1.d-4

end module amrex_eb_flux_reg_nd_module
