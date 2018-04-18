
module amrex_bc_types_module

  implicit none

!  include 'AMReX_bc_types.fi'
  
  private

  integer, parameter, public :: amrex_bc_bogus        = -666
  integer, parameter, public :: amrex_bc_reflect_odd  = -1
  integer, parameter, public :: amrex_bc_int_dir      =  0
  integer, parameter, public :: amrex_bc_reflect_even =  1
  integer, parameter, public :: amrex_bc_foextrap     =  2
  integer, parameter, public :: amrex_bc_ext_dir      =  3
  integer, parameter, public :: amrex_bc_hoextrap     =  4

end module amrex_bc_types_module
