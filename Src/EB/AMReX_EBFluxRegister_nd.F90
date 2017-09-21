module amrex_eb_flux_reg_nd_module
  use amrex_fort_module, only : amrex_real
  implicit none
  private

  integer, parameter, public :: crse_cell = 0
  integer, parameter, public :: crse_fine_boundary_cell = 1
  integer, parameter, public :: fine_cell = 2
  real(amrex_real), parameter, public :: reredistribution_threshold = 1.d-4

end module amrex_eb_flux_reg_nd_module
