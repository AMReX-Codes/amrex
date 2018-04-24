
!-----------------------------------------------------------------------

  subroutine bl_proffortfuncstart(str)

    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)

    call bl_str2int(istr, NSTR, str)
    call bl_proffortfuncstart_cpp(istr, NSTR)

  end subroutine bl_proffortfuncstart

!-----------------------------------------------------------------------

  subroutine bl_proffortfuncstop(str)

    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)

    call bl_str2int(istr, NSTR, str)
    call bl_proffortfuncstop_cpp(istr, NSTR)

  end subroutine bl_proffortfuncstop

!-----------------------------------------------------------------------

  subroutine bl_proffortfuncstart_int(i)

    integer :: i

    call bl_proffortfuncstart_cpp_int(i)

  end subroutine bl_proffortfuncstart_int

!-----------------------------------------------------------------------

  subroutine bl_proffortfuncstop_int(i)

    integer :: i

    call bl_proffortfuncstop_cpp_int(i)

  end subroutine bl_proffortfuncstop_int
