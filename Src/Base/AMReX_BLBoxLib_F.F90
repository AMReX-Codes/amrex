
!-----------------------------------------------------------------------

  subroutine bl_error(str)

    character*(*) :: str
    integer, parameter :: NSTR = 256
    integer :: istr(NSTR)

    flush(6)
    call blstr2int(istr, NSTR, str)
    call bl_error_cpp(istr, NSTR)

  end subroutine bl_error

!-----------------------------------------------------------------------

  subroutine bl_warning(str)

    character*(*) :: str
    integer, parameter :: NSTR = 256
    integer :: istr(NSTR)

    flush(6)
    call blstr2int(istr, NSTR, str)
    call bl_warning_cpp(istr, NSTR)

  end subroutine bl_warning

!-----------------------------------------------------------------------

  subroutine bl_abort(str)

    character*(*) :: str
    integer, parameter :: NSTR = 256
    integer :: istr(NSTR)

    flush(6)
    call blstr2int(istr, NSTR, str)
    call bl_abort_cpp(istr, NSTR)

  end subroutine bl_abort
