
module amrex_bc_types_module

  implicit none

  include 'AMReX_bc_types.fi'
  
  private

  integer, parameter, public :: amrex_bc_bogus        = BOGUS_BC
  integer, parameter, public :: amrex_bc_reflect_odd  = REFLECT_ODD 
  integer, parameter, public :: amrex_bc_int_dir      = INT_DIR
  integer, parameter, public :: amrex_bc_reflect_even = REFLECT_EVEN
  integer, parameter, public :: amrex_bc_foextrap     = FOEXTRAP
  integer, parameter, public :: amrex_bc_ext_dir      = EXT_DIR
  integer, parameter, public :: amrex_bc_hoextrap     = HOEXTRAP

end module amrex_bc_types_module
