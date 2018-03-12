
module warpx_fft_module
  use amrex_fort_module
  implicit none

  private
  public :: warpx_fft_mpi_init

contains

  subroutine warpx_fft_mpi_init (fcomm) bind(c,name='warpx_fft_mpi_init')
    use shared_data, only : comm, rank, nproc
    integer, intent(in), value :: fcomm

    integer :: ierr, lnproc, lrank
  
    comm = fcomm

    call mpi_comm_size(comm, lnproc, ierr)
    nproc = lnproc

    call mpi_comm_rank(comm, lrank, ierr)
    rank = lrank
  end subroutine warpx_fft_mpi_init

end module warpx_fft_module
