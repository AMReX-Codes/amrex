
subroutine bar (comm) bind(c)
  use amrex_base_module
  integer, intent(in), value :: comm

  integer :: rank, ierr

  call amrex_init(comm)

  if (amrex_parallel_ioprocessor()) then
     print *, "bar: AMReX Fortran has been initialized."
  end if

  call amrex_finalize()

  ! After amrex_finalize(), amrex can no longer be used.
  call mpi_comm_rank(comm, rank, ierr)
  if (rank .eq. 0) then
     print *, "bar: AMReX Fortran has been finalized."
  end if
  
end subroutine bar
