
!-----------------------------------------------------------------------

  subroutine bl_pp_new(ipp, str)

    integer :: ipp
    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)

    call blstr2int(istr, NSTR, str)
    call bl_pp_new_cpp(ipp, istr, NSTR)

  end subroutine bl_pp_new

!-----------------------------------------------------------------------

  subroutine bl_pp_record_new(ipp, ippr, str)

    integer :: ipp, ippr
    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)

    call blstr2int(istr, NSTR, str)
    call bl_pp_record_new_cpp(ipp, ippr, istr, NSTR)

  end subroutine bl_pp_record_new

!-----------------------------------------------------------------------

  subroutine bl_pp_get_int(ierr, ipp, str, ival)

    integer :: ierr, ipp
    integer :: ival
    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)
    
    call blstr2int(istr, NSTR, str)
    call bl_pp_get_int_cpp(ierr, ipp, istr, NSTR, ival)

  end subroutine bl_pp_get_int

!-----------------------------------------------------------------------

  subroutine bl_pp_get_int_n(ierr, ipp, str, ival, nval)

    integer :: ierr, ipp, nval
    integer :: ival(*)
    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)

    call blstr2int(istr, NSTR, str)
    call bl_pp_get_int_n_cpp(ierr, ipp, istr, NSTR, ival, nval)

  end subroutine bl_pp_get_int_n

!-----------------------------------------------------------------------

  subroutine bl_pp_get_logical(ierr, ipp, str, lval)

    integer :: ierr, ipp
    logical :: lval
    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)
    integer :: ival

    call blstr2int(istr, NSTR, str)
    call bl_pp_get_logical_cpp(ierr, ipp, istr, NSTR, ival)

    if ( IERR .NE. 0 ) then
       if ( ival .NE. 0 ) then
          lval = .TRUE.
       else
          lval = .FALSE.
       end if
    end if

  end subroutine bl_pp_get_logical

!-----------------------------------------------------------------------

  subroutine bl_pp_get_real(ierr, ipp, str, rval)

    integer :: ierr, ipp
    real :: rval
    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)

    call blstr2int(istr, NSTR, str)
    call bl_pp_get_real_cpp(ierr, ipp, istr, NSTR, rval)

  end subroutine bl_pp_get_real

!-----------------------------------------------------------------------

  subroutine bl_pp_get_real_n(ierr, ipp, str, rval,nval)

    integer :: ierr, ipp, nval
    real :: rval(*)
    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)

    call blstr2int(istr, NSTR, str)
    call bl_pp_get_real_n_cpp(ierr, ipp, istr, NSTR, rval,nval)

  end subroutine bl_pp_get_real_n

!-----------------------------------------------------------------------

  subroutine bl_pp_get_double(ierr, ipp, str, dval)

    integer :: ierr, ipp
    double precision :: dval
    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)

    call blstr2int(istr, NSTR, str)
    call bl_pp_get_double_cpp(ierr, ipp, istr, NSTR, dval)

  end subroutine bl_pp_get_double

!-----------------------------------------------------------------------

  subroutine bl_pp_get_double_n(ierr, ipp, str, dval, nval)

    integer :: ierr, ipp, nval
    double precision :: dval
    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)

    call blstr2int(istr, NSTR, str)
    call bl_pp_get_double_n_cpp(ierr, ipp, istr, NSTR, dval,nval)

  end subroutine bl_pp_get_double_n

!-----------------------------------------------------------------------

  subroutine bl_pp_get_string(ierr, ipp, str, ostr)

    integer :: ierr, ipp
    character*(*) :: ostr
    character*(*) :: str
    integer, parameter :: NSTR = 128
    integer :: istr(NSTR)
    integer :: iostr(NSTR)

    call blstr2int(istr, NSTR, str)
    call bl_pp_get_string_cpp(ierr, ipp, istr, NSTR, iostr, NSTR)

    if ( ierr .ne. 0 ) then
       call blint2str(ostr, iostr, NSTR)
    end if

  end subroutine bl_pp_get_string
