!!
!! Provide methods for accessing problem definition files and the command line
!!
module bl_parmparse_module

  implicit none

contains

  !! Seeks to a position in an input file corresponding to the
  !! start of an NAMELIST record.
  !! The '&NML' must be a the beginning of the line, NML data
  !! items can follow the &NML on the line, however.
  !! On exit, the I/O UNIT is positioned at the so that a read
  !! on that namelist will succeed
  !! Else the I/O UNIT will be at EOF
  !! If STAT is not present, the BL_ERROR will be called on EOF, and
  !! the NML is not found.

  subroutine parmparse_seek_to(unit, nml, stat)

    use bl_error_module

    integer, intent(in) :: unit
    character(len=*), intent(in) :: nml
    logical, intent(out), optional :: stat

    integer, parameter :: MAX_INPUTS_LEN   = 256
    character(len=len(nml)+1)    :: cnml
    character(len=MAX_INPUTS_LEN) :: buf
    integer :: ierr, cnml_len

    cnml = "&" // nml
    cnml_len = len_trim(cnml)
    rewind(unit = unit)
    do
       read(unit=unit, fmt='(a)', end = 100, iostat = ierr) buf
       if ( ierr > 0 ) then
          call bl_error("PARMPARSE_SEEK_TO: IO Error on parse inputs unit ", unit)
       end if
       if ( cnml == buf(1:cnml_len) ) then
          backspace(unit=unit)
          if ( present(stat) ) stat = .TRUE.
          return
       end if
    end do
100 continue
    if ( present(stat) ) then
       stat = .FALSE.
       return
    else
       call bl_error("PARMPARSE_SEEK_TO: not found namelist : ", nml)
    end if

  end subroutine parmparse_seek_to

  !! Returns the number of arguments on the command line that match
  !! STR, or if present STR_LONG.
  function pp_arg_count(str, str_long) result(r)
    use f2kcli
    character(len=*), intent(in) :: str
    character(len=*), intent(in), optional :: str_long
    integer :: r
    integer :: narg, f
    character(len=128) :: option
    narg = command_argument_count()
    r = 0
    do f = 1, narg
       call get_command_argument(f, value = option )
       if ( str == option ) then
          r = r + 1
       else if ( present(str_long) ) then
          if ( str_long == option ) then
             r = r + 1
          end if
       end if
    end do
  end function pp_arg_count

end module bl_parmparse_module
