module amrex_ya_flux_reg_nd_module
  use amrex_fort_module, only : amrex_real
  implicit none
  private

  integer, parameter, public :: crse_cell = 0
  integer, parameter, public :: crse_fine_boundary_cell = 1
  integer, parameter, public :: fine_cell = 2

end module amrex_ya_flux_reg_nd_module
