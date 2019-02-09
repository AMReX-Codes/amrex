module amrex_IO_module

  implicit none

  !> Default input and output units:

  integer, private, parameter :: IO_STDIN  = 5
  integer, private, parameter :: IO_STDOUT = 6

  !> Number and value of pre-connected units

  integer, private, parameter :: IO_NUM_PRECON = 3
  integer, private, parameter :: IO_PRECON_UNITS(IO_NUM_PRECON) = (/ 0, 5, 6 /)

  !> Largest allowed unit number (or a large number, if none)

  integer, private, parameter :: IO_MAX_UNIT = 1000

  ! Note that the constants defined in this module are not portable,
  ! and may have to be changed for any given processor.

  !> Several routines are written so that they can be conveniently
  !! called in routines that have optional arguments.  For example,
  !! ``
  !! subroutine print_box(box, unit)
  !!  optional :: unit
  !!  write(unit=unit_stdout(unit), fmt *) box
  !! ``
  !! This routine will print on the stdout if UNIT is not passed to
  !! print_box, but will print to UNIT if UNIT is passed.

contains

  !> Returns a unit number of a unit that exists and is not connected
  function unit_new() result (r)
    integer :: r
    logical :: exists, opened
    integer :: ios
    do r = 1, IO_MAX_UNIT
       if (r == IO_STDIN .or. r == IO_STDOUT ) cycle
       if (any (r == IO_PRECON_UNITS )) cycle
       inquire (unit = r,  &
            exist = exists,  &
            opened = opened,  &
            iostat = ios)
       if ( exists .and. .not. opened .and. ios == 0 ) return
    end do
    r = -1
  end function unit_new

  !> Returns the stdin unit number if no argument is passed, or else
  !! unit if it is passed.
  pure function unit_stdin(unit)  result (r)
    integer :: r
    integer, intent(in), optional :: unit
    if ( present(unit) ) then
       r = unit
    else
       r = IO_STDIN
    end if
  end function unit_stdin

  !> Returns the stdout unit number if no argument is passed, or else
  !! unit if it is passed.
  pure function unit_stdout(unit) result(r)
    integer :: r
    integer, intent(in), optional :: unit
    if ( present(unit) ) then
       r = unit
    else
       r = IO_STDOUT
    end if
  end function unit_stdout

  !> Returns the string 'YES' if no argument is passed other wise
  !! returns the argument advance
  pure function unit_advance(advance) result(r)
    character(len=3) :: r
    character(len=*), intent(in), optional :: advance
    if ( present(advance) ) then
       r = advance
    else
       r = 'YES'
    end if
  end function unit_advance

  !> Puts skip spaces of output, without advancing to the next
  !! record onto UNIT.  If no skip argument is given, then no advance
  !! is done.
  subroutine unit_skip(unit, skip)
    integer, intent(in) :: unit
    integer, intent(in), optional :: skip
    if ( .not. present(skip) ) return
    write(unit=unit, fmt='(A)', advance = 'NO') repeat(' ', skip)
  end subroutine unit_skip

  !> A convenience function that returns 0 if SKIP is not present, otherwise
  !! it returns SKIP
  pure function unit_get_skip(skip) result(r)
    integer :: r
    integer, intent(in), optional :: skip
    r = 0; if ( present(skip) ) r = skip
  end function unit_get_skip

end module amrex_IO_module
