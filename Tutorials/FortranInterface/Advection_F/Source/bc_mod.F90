
module bc_module

  use amrex_base_module

  implicit none

  ! periodic bc.  See amrex_bc_types_module for a list of bc types.
  integer, save :: lo_bc(amrex_spacedim,1) = amrex_bc_int_dir  ! the second dimension is the
  integer, save :: hi_bc(amrex_spacedim,1) = amrex_bc_int_dir  ! number of components
  
end module bc_module
