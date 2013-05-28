!! Support routines for the _BoxLib_ framework.
!!
!! The _BoxLib_ framework provides access to the parallel
!! computing environment, error reporting routings, defines
!! some basic types and parameter.
module BoxLib

  use parallel
  use omp_module

  implicit none

  private

  public :: boxlib_initialize
  public :: boxlib_finalize

contains

  !! Initializes _BoxLib_ applications.  This should be the
  !! first routine called in the main PROGRAM unit.
  subroutine boxlib_initialize(thread_support_level)
    integer, intent(in), optional :: thread_support_level
    call parallel_initialize(MPI_COMM_WORLD, thread_support_level)
    if (parallel_IOProcessor() .and. parallel_nprocs() > 1) then
       print*, "MPI initialized with ", parallel_nprocs(), " MPI processes";
       if (present(thread_support_level)) then
          select case(parallel_thread_support_level())
          case (MPI_THREAD_SINGLE)
             print*, "MPI thread support level: MPI_THREAD_SINGLE"
          case (MPI_THREAD_FUNNELED)
             print*, "MPI thread support level: MPI_THREAD_FUNNELED"
          case (MPI_THREAD_SERIALIZED)
             print*, "MPI thread support level: MPI_THREAD_SERIALIZED"
          case (MPI_THREAD_MULTIPLE)
             print*, "MPI thread support level: MPI_THREAD_MULTIPLE"
          case default
             print*, "WARNING: Unknown MPI thread support level: ", parallel_thread_support_level()
          end select
       end if
       print*, "MPI initialized with ", omp_get_max_threads(), " threads";
    endif
    if (omp_get_max_threads() > 1) call omp_set_nested(.false.)
  end subroutine boxlib_initialize

  !! Finalizes _BoxLib_ applications. This should be the final
  !! routine called in the main PROGRAM unit.
  subroutine boxlib_finalize()
    call parallel_finalize
  end subroutine boxlib_finalize

end module BoxLib
