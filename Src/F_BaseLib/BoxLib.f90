!! Support routines for the _BoxLib_ framework.
!!
!! The _BoxLib_ framework provides access to the parallel
!! computing environment, error reporting routings, defines
!! some basic types and parameter.
module BoxLib

  use parallel

  implicit none

  private

  public :: boxlib_initialize
  public :: boxlib_finalize

contains

  !! Initializes _BoxLib_ applications.  This should be the
  !! first routine called in the main PROGRAM unit.
  subroutine boxlib_initialize()
    call parallel_initialize()
  end subroutine boxlib_initialize

  !! Finalizes _BoxLib_ applications. This should be the final
  !! routine called in the main PROGRAM unit.
  subroutine boxlib_finalize()
    call parallel_finalize
  end subroutine boxlib_finalize

end module BoxLib
