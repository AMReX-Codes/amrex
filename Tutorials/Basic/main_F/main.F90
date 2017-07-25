
program main

  use mpi
  use amrex_base_module

  implicit none

  integer :: ierr

  call mpi_init(ierr)

  call amrex_init()

  
  

  call amrex_finalize()

  call mpi_finalize()

end program main
