
module bc_module

  use amrex_base_module

  implicit none

  ! periodic bc.  See amrex_bc_types_module for a list of bc types.
  integer, parameter :: lo_bc(amrex_spacedim,1) = amrex_bc_int_dir  ! the second dimension is the
  integer, parameter :: hi_bc(amrex_spacedim,1) = amrex_bc_int_dir  ! number of components
  
contains

end module bc_module
