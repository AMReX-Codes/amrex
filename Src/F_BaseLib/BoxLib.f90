!! Support routines for the _BoxLib_ framework.
!!
!! The _BoxLib_ framework provides access to the parallel
!! computing environment, error reporting routings, defines
!! some basic types and parameter.
module BoxLib

  use bl_error_module
  use bl_space
  use bl_types
  use parallel

  implicit none

  !! The default dump unint
  integer, private :: dunit = -1

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

  !! Returns the dump unit.
  function boxlib_dunit() result(r)
    integer :: r
    r = dunit
  end function boxlib_dunit

  !! Sets the dump unit.
  subroutine boxlib_set_dunit(unit)
    integer, intent(in) :: unit
    dunit = unit
  end subroutine boxlib_set_dunit

  !! Sets the Dump outfile.  Without setting this file
  subroutine boxlib_open_dfile(str)
    use bl_IO_module
    character(len=*) str
    dunit = unit_new()
    open(unit=dunit, file = str, status='replace', action = 'write')
  end subroutine boxlib_open_dfile

  !! Closes the dump unit.  Will no longer be available for IO.
  subroutine boxlib_close_dfile()
    use bl_IO_module
    close(unit=dunit)
    dunit = -1
  end subroutine boxlib_close_dfile

end module BoxLib
