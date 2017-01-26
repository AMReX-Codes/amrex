
#ifdef BL_USE_MPI
module fboxlib_mpi
  implicit none
  include 'mpif.h'
end module fboxlib_mpi
#endif

module parallel_module

  use iso_c_binding
#ifdef BL_USE_MPI
  use fboxlib_mpi
#endif

  implicit none

  private

  public :: parallel_comm_init_from_c
  public :: parallel_comm_free_from_c
  public :: parallel_myproc 
  public :: parallel_ioprocessor
  public :: parallel_nprocs

#ifdef BL_USE_MPI
  integer :: m_nprocs = -1
  integer :: m_myproc = -1
  integer :: m_comm   = -1
#else
  integer :: m_nprocs = 0
  integer :: m_myproc = 0
  integer :: m_comm   = 0
#endif

contains

  subroutine parallel_comm_init_from_c (comm) bind(c, name='bl_fortran_mpi_comm_init')
    use iso_c_binding
    integer(c_int), intent(in), value :: comm
#ifdef BL_USE_MPI
    integer :: ierr
    call MPI_Comm_Dup(comm, m_comm, ierr)
    call MPI_Comm_Size(m_comm, m_nprocs, ierr)
    call MPI_Comm_Rank(m_comm, m_myproc, ierr)
    call MPI_barrier(m_comm, ierr)
#else
    m_nprocs = 1
    m_myproc = 0
#endif
  end subroutine parallel_comm_init_from_c

  subroutine parallel_comm_free_from_c () bind(c, name='bl_fortran_mpi_comm_free')
#ifdef BL_USE_MPI
    integer :: ierr
    call MPI_Comm_Free(m_comm, ierr)
#endif
  end subroutine parallel_comm_free_from_c

  integer function parallel_myproc ()
    parallel_myproc = m_myproc
  end function parallel_myproc

  logical function parallel_ioprocessor ()
    parallel_ioprocessor = m_myproc .eq. 0
  end function parallel_ioprocessor

  integer function parallel_nprocs ()
    parallel_nprocs = m_nprocs
  end function parallel_nprocs

end module parallel_module
