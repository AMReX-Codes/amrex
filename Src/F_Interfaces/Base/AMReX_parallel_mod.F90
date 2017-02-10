
#ifdef BL_USE_MPI
module amrex_fi_mpi
  implicit none
  include 'mpif.h'
end module amrex_fi_mpi
#endif

module amrex_parallel_module

  use iso_c_binding
#ifdef BL_USE_MPI
  use amrex_fi_mpi
#endif

  implicit none

  private

  public :: amrex_parallel_comm_init_from_c
  public :: amrex_parallel_comm_free_from_c
  public :: amrex_parallel_myproc 
  public :: amrex_parallel_ioprocessor
  public :: amrex_parallel_nprocs

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

  subroutine amrex_parallel_comm_init_from_c (comm) bind(c, name='amrex_parallel_comm_init_from_c')
    use iso_c_binding
    integer(c_int), value :: comm
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
  end subroutine amrex_parallel_comm_init_from_c

  subroutine amrex_parallel_comm_free_from_c () bind(c, name='amrex_parallel_comm_free_from_c')
#ifdef BL_USE_MPI
    integer :: ierr
    call MPI_Comm_Free(m_comm, ierr)
#endif
  end subroutine amrex_parallel_comm_free_from_c

  integer function amrex_parallel_myproc ()
    amrex_parallel_myproc = m_myproc
  end function amrex_parallel_myproc

  logical function amrex_parallel_ioprocessor ()
    amrex_parallel_ioprocessor = m_myproc .eq. 0
  end function amrex_parallel_ioprocessor

  integer function amrex_parallel_nprocs ()
    amrex_parallel_nprocs = m_nprocs
  end function amrex_parallel_nprocs

end module amrex_parallel_module
