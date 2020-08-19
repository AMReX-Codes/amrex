
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

  use amrex_error_module
  use amrex_fort_module, only : amrex_real

  implicit none

  private

  public :: amrex_parallel_init
  public :: amrex_parallel_finalize
  public :: amrex_parallel_communicator
  public :: amrex_parallel_myproc 
  public :: amrex_parallel_ioprocessor
  public :: amrex_parallel_nprocs
  public :: amrex_parallel_reduce_sum
  public :: amrex_parallel_reduce_max
  public :: amrex_parallel_reduce_min

#ifdef BL_USE_MPI
  integer, public :: amrex_mpi_real = MPI_DATATYPE_NULL
  integer :: m_nprocs = -1
  integer :: m_myproc = -1
  integer :: m_comm   = MPI_COMM_NULL
#else
  integer :: m_nprocs = 0
  integer :: m_myproc = 0
  integer :: m_comm   = 0
#endif

  interface amrex_parallel_reduce_sum
     module procedure amrex_parallel_reduce_sum_is
     module procedure amrex_parallel_reduce_sum_iv
     module procedure amrex_parallel_reduce_sum_rs
     module procedure amrex_parallel_reduce_sum_rv
  end interface amrex_parallel_reduce_sum

  interface amrex_parallel_reduce_max
     module procedure amrex_parallel_reduce_max_is
     module procedure amrex_parallel_reduce_max_iv
     module procedure amrex_parallel_reduce_max_rs
     module procedure amrex_parallel_reduce_max_rv
  end interface amrex_parallel_reduce_max

  interface amrex_parallel_reduce_min
     module procedure amrex_parallel_reduce_min_is
     module procedure amrex_parallel_reduce_min_iv
     module procedure amrex_parallel_reduce_min_rs
     module procedure amrex_parallel_reduce_min_rv
  end interface amrex_parallel_reduce_min

  logical, save :: call_mpi_finalize = .false.

contains

  subroutine amrex_parallel_init (comm)
    integer, intent(in), optional :: comm
#ifdef BL_USE_MPI
    integer :: ierr
    logical :: flag
    call MPI_Initialized(flag, ierr)

    if (present(comm) .and. .not.flag) then
       if (comm .ne. MPI_COMM_WORLD) then
          stop "MPI has not been initialized.  How come we are given a communciator?"
       endif
    end if

    if (.not.flag) then
       call MPI_Init(ierr)
       call_mpi_finalize = .true.
    else
       call_mpi_finalize = .false.
    end if

    if (present(comm)) then
       call MPI_Comm_Dup(comm, m_comm, ierr)
    else
       call MPI_Comm_Dup(MPI_COMM_WORLD, m_comm, ierr)
    end if

    call MPI_Comm_Size(m_comm, m_nprocs, ierr)
    call MPI_Comm_Rank(m_comm, m_myproc, ierr)
    if (amrex_real == c_double) then
       amrex_mpi_real = MPI_DOUBLE_PRECISION
    else if (amrex_real == c_float) then
       amrex_mpi_real = MPI_REAL
    else
       call amrex_abort("amrex_parallel_init: size of amrex_real is unknown")
    end if
#endif    
  end subroutine amrex_parallel_init

  subroutine amrex_parallel_finalize ()
#ifdef BL_USE_MPI
    integer :: ierr
    call MPI_Comm_Free(m_comm, ierr)
    m_comm = MPI_COMM_NULL
    if (call_mpi_finalize) then
       call MPI_Finalize(ierr)
       call_mpi_finalize = .false.
    end if
#else
    call_mpi_finalize = .false.
#endif
  end subroutine amrex_parallel_finalize

  pure integer function amrex_parallel_communicator ()
    amrex_parallel_communicator = m_comm
  end function amrex_parallel_communicator

  pure integer function amrex_parallel_myproc ()
    amrex_parallel_myproc = m_myproc
  end function amrex_parallel_myproc

  pure logical function amrex_parallel_ioprocessor ()
    amrex_parallel_ioprocessor = m_myproc .eq. 0
  end function amrex_parallel_ioprocessor

  pure integer function amrex_parallel_nprocs ()
    amrex_parallel_nprocs = m_nprocs
  end function amrex_parallel_nprocs

  subroutine amrex_parallel_reduce_sum_is (i, rank)
    integer, intent(inout) :: i
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    integer :: tmp, ierr
    tmp = i
    if (present(rank)) then
       call MPI_Reduce(tmp, i, 1, MPI_INTEGER, MPI_SUM, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, i, 1, MPI_INTEGER, MPI_SUM, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_sum_is

  subroutine amrex_parallel_reduce_sum_iv (i, n, rank)
    integer, intent(inout) :: i(*)
    integer, intent(in) :: n
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    integer :: tmp(n), ierr
    tmp = i(1:n)
    if (present(rank)) then
       call MPI_Reduce(tmp, i, n, MPI_INTEGER, MPI_SUM, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, i, n, MPI_INTEGER, MPI_SUM, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_sum_iv

  subroutine amrex_parallel_reduce_sum_rs (r, rank)
    real(amrex_real), intent(inout) :: r
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    real(amrex_real) :: tmp
    integer :: ierr
    tmp = r
    if (present(rank)) then
       call MPI_Reduce(tmp, r, 1, amrex_mpi_real, MPI_SUM, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, r, 1, amrex_mpi_real, MPI_SUM, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_sum_rs

  subroutine amrex_parallel_reduce_sum_rv (r, n, rank)
    real(amrex_real), intent(inout) :: r(*)
    integer, intent(in) :: n
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    real(amrex_real) :: tmp(n)
    integer :: ierr
    tmp = r(1:n)
    if (present(rank)) then
       call MPI_Reduce(tmp, r, n, amrex_mpi_real, MPI_SUM, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, r, n, amrex_mpi_real, MPI_SUM, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_sum_rv

  subroutine amrex_parallel_reduce_max_is (i, rank)
    integer, intent(inout) :: i
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    integer :: tmp, ierr
    tmp = i
    if (present(rank)) then
       call MPI_Reduce(tmp, i, 1, MPI_INTEGER, MPI_MAX, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, i, 1, MPI_INTEGER, MPI_MAX, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_max_is

  subroutine amrex_parallel_reduce_max_iv (i, n, rank)
    integer, intent(inout) :: i(*)
    integer, intent(in) :: n
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    integer :: tmp(n), ierr
    tmp = i(1:n)
    if (present(rank)) then
       call MPI_Reduce(tmp, i, n, MPI_INTEGER, MPI_MAX, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, i, n, MPI_INTEGER, MPI_MAX, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_max_iv

  subroutine amrex_parallel_reduce_max_rs (r, rank)
    real(amrex_real), intent(inout) :: r
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    real(amrex_real) :: tmp
    integer :: ierr
    tmp = r
    if (present(rank)) then
       call MPI_Reduce(tmp, r, 1, amrex_mpi_real, MPI_MAX, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, r, 1, amrex_mpi_real, MPI_MAX, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_max_rs

  subroutine amrex_parallel_reduce_max_rv (r, n, rank)
    real(amrex_real), intent(inout) :: r(*)
    integer, intent(in) :: n
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    real(amrex_real) :: tmp(n)
    integer :: ierr
    tmp = r(1:n)
    if (present(rank)) then
       call MPI_Reduce(tmp, r, n, amrex_mpi_real, MPI_MAX, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, r, n, amrex_mpi_real, MPI_MAX, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_max_rv

  subroutine amrex_parallel_reduce_min_is (i, rank)
    integer, intent(inout) :: i
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    integer :: tmp, ierr
    tmp = i
    if (present(rank)) then
       call MPI_Reduce(tmp, i, 1, MPI_INTEGER, MPI_MIN, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, i, 1, MPI_INTEGER, MPI_MIN, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_min_is

  subroutine amrex_parallel_reduce_min_iv (i, n, rank)
    integer, intent(inout) :: i(*)
    integer, intent(in) :: n
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    integer :: tmp(n), ierr
    tmp = i(1:n)
    if (present(rank)) then
       call MPI_Reduce(tmp, i, n, MPI_INTEGER, MPI_MIN, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, i, n, MPI_INTEGER, MPI_MIN, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_min_iv

  subroutine amrex_parallel_reduce_min_rs (r, rank)
    real(amrex_real), intent(inout) :: r
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    real(amrex_real) :: tmp
    integer :: ierr
    tmp = r
    if (present(rank)) then
       call MPI_Reduce(tmp, r, 1, amrex_mpi_real, MPI_MIN, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, r, 1, amrex_mpi_real, MPI_MIN, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_min_rs

  subroutine amrex_parallel_reduce_min_rv (r, n, rank)
    real(amrex_real), intent(inout) :: r(*)
    integer, intent(in) :: n
    integer, intent(in), optional :: rank
#ifdef BL_USE_MPI
    real(amrex_real) :: tmp(n)
    integer :: ierr
    tmp = r(1:n)
    if (present(rank)) then
       call MPI_Reduce(tmp, r, n, amrex_mpi_real, MPI_MIN, rank, m_comm, ierr)
    else
       call MPI_Allreduce(tmp, r, n, amrex_mpi_real, MPI_MIN, m_comm, ierr)
    end if
#endif
  end subroutine amrex_parallel_reduce_min_rv

end module amrex_parallel_module
