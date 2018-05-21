
module bc_module

  use amrex_base_module

  implicit none

  ! periodic bc.  See amrex_bc_types_module for a list of bc types.
  integer, allocatable :: lo_bc(:,:)
  integer, allocatable :: hi_bc(:,:)
  
end module bc_module
