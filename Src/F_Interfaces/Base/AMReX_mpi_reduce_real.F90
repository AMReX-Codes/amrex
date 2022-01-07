! Starting from gfortran 10, mismatches between actual and dummy argument
! lists in a single file have been rejected with an error.  This causes
! issues for mpich. https://lists.mpich.org/pipermail/discuss/2020-January/005863.html
! This is a workaround by splitting calls to int and real into two files.

module amrex_mpi_reduce_real_module
  use amrex_fi_mpi
  use amrex_fort_module, only : amrex_real
  implicit none
  private
  public :: amrex_mpi_reduce_real, amrex_mpi_allreduce_real

  interface amrex_mpi_reduce_real
     module procedure amrex_mpi_reduce_real_s
     module procedure amrex_mpi_reduce_real_v
  end interface amrex_mpi_reduce_real

  interface amrex_mpi_allreduce_real
     module procedure amrex_mpi_allreduce_real_s
     module procedure amrex_mpi_allreduce_real_v
  end interface amrex_mpi_allreduce_real

contains

  subroutine amrex_mpi_reduce_real_s (sendbuf, recvbuf, count, datatype, op, root, comm, ierror)
    real(amrex_real), intent(in) :: sendbuf
    real(amrex_real), intent(out) :: recvbuf
    integer, intent(in) :: count, datatype, op, root, comm
    integer, intent(out) :: ierror
    real(amrex_real) :: src(1), dst(1)
    src(1) = sendbuf
    call MPI_Reduce(src, dst, 1, datatype, op, root, comm, ierror)
    recvbuf = dst(1)
  end subroutine amrex_mpi_reduce_real_s

  subroutine amrex_mpi_reduce_real_v (sendbuf, recvbuf, count, datatype, op, root, comm, ierror)
    real(amrex_real), intent(in) :: sendbuf(*)
    real(amrex_real) :: recvbuf(*)
    integer, intent(in) :: count, datatype, op, root, comm
    integer, intent(out) :: ierror
    call MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm, ierror)
  end subroutine amrex_mpi_reduce_real_v

  subroutine amrex_mpi_allreduce_real_s (sendbuf, recvbuf, count, datatype, op, comm, ierror)
    real(amrex_real), intent(in) :: sendbuf
    real(amrex_real), intent(out) :: recvbuf
    integer, intent(in) :: count, datatype, op, comm
    integer, intent(out) :: ierror
    real(amrex_real) :: src(1), dst(1)
    src(1) = sendbuf
    call MPI_Allreduce(src, dst, 1, datatype, op, comm, ierror)
    recvbuf = dst(1)
  end subroutine amrex_mpi_allreduce_real_s

  subroutine amrex_mpi_allreduce_real_v (sendbuf, recvbuf, count, datatype, op, comm, ierror)
    real(amrex_real), intent(in) :: sendbuf(*)
    real(amrex_real) :: recvbuf(*)
    integer, intent(in) :: count, datatype, op, comm
    integer, intent(out) :: ierror
    call MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierror)
  end subroutine amrex_mpi_allreduce_real_v

end module amrex_mpi_reduce_real_module
