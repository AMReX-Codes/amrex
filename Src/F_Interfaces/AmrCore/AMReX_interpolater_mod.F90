
module amrex_interpolater_module
  implicit none
  ! THIS MUST BE CONSISTENT WITH amrex_fi_fillpatch_two in AMReX_fillpatch_fi.cpp!!!
  integer, parameter :: amrex_interp_pc            = 0
  integer, parameter :: amrex_interp_node_bilinear = 1
  integer, parameter :: amrex_interp_cell_bilinear = 2
  integer, parameter :: amrex_interp_quadratic     = 3
  integer, parameter :: amrex_interp_lincc         = 4
  integer, parameter :: amrex_interp_cell_cons     = 5
  integer, parameter :: amrex_interp_protected     = 6
  integer, parameter :: amrex_interp_quartic       = 7
end module amrex_interpolater_module
