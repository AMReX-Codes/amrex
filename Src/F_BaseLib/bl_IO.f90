module bl_IO_module

  implicit none

  ! Values returned to IOSTAT for end of record and end of file

  integer, parameter :: END_OF_RECORD = -2
  integer, parameter :: END_OF_FILE = -1

  ! Default input and output units:

  integer, parameter :: DEFAULT_INPUT_UNIT = 5
  integer, parameter :: DEFAULT_OUTPUT_UNIT = 6

  ! Number and value of pre-connected units

  integer, parameter :: NUMBER_OF_PRECONNECTED_UNITS = 3
  integer, parameter :: PRECONNECTED_UNITS (NUMBER_OF_PRECONNECTED_UNITS) = &
       (/ 0, 5, 6 /)

  ! Largest allowed unit number (or a large number, if none)

  integer, parameter :: MAX_UNIT_NUMBER = 1000

  integer, parameter :: MAX_INPUT_LEN   = 1024

contains

  function unit_new ()  result (r)

    ! Returns a unit number of a unit that exists and is not connected

    integer :: r
    logical :: exists, opened
    integer :: ios

    do r = 1, max_unit_number
       if (r == DEFAULT_INPUT_UNIT .or. &
            r == DEFAULT_OUTPUT_UNIT) cycle
       if (any (r == PRECONNECTED_UNITS)) cycle
       inquire (unit = r,  &
            exist = exists,  &
            opened = opened,  &
            iostat = ios)
       if (exists .and. .not. opened .and. ios == 0) return
    end do

    r = -1

  end function unit_new

  function unit_stdin (unit)  result (r)
    integer :: r
    integer, intent(in), optional :: unit
    if (present(unit)) then
       r = unit
    else
       r = DEFAULT_INPUT_UNIT
    end if
  end function unit_stdin

  function unit_stdout (unit) result(r)
    integer :: r
    integer, intent(in), optional :: unit
    if ( present(unit) ) then
       r = unit
    else
       r = DEFAULT_OUTPUT_UNIT
    end if
  end function unit_stdout

  function unit_advance (advance) result(r)
    character(len=3) :: r
    character(len=*), intent(in), optional :: advance
    if ( present(advance) ) then
       r = advance
    else
       r = 'YES'
    end if
  end function unit_advance

  subroutine unit_skip(unit, skip)
    integer, intent(in) :: unit
    integer, intent(in), optional :: skip
    integer :: i
    if ( .not. present(skip) ) return
    do i = 1, skip
       write(unit=unit, fmt='(" ")', advance = 'NO')
    end do
  end subroutine unit_skip

  function unit_get_skip(skip) result(r)
    integer :: r
    integer, intent(in), optional :: skip
    r = 0; if ( present(skip) ) r = skip
  end function unit_get_skip

  subroutine inputs_seek_to(unit, nml, stat)
    use bl_error_module
    character(len=*), intent(in) :: nml
    integer, intent(in) :: unit
    integer, intent(out), optional :: stat
    character(len=MAX_INPUT_LEN) :: buf
    character(len=len(nml)+1)    :: cnml
    integer :: ierr, cnml_len

    cnml = "&" // nml
    cnml_len = len_trim(cnml)
    rewind(unit = unit)
    do
       read(unit=unit, fmt='(a)', end = 100, iostat = ierr) buf
       if ( ierr > 0 ) then
          call bl_error("INPUTS_SEEK_TO: IO Error on inputs: stat = ", ierr)
       end if
       if ( cnml == buf(1:cnml_len) ) then
          backspace(unit=unit)
          if ( present(stat) ) stat = 0
          return
       end if
    end do
100 continue
    if ( present(stat) ) then
       stat = 1
       return
    end if
    call bl_error("INPUTS_SEEK_TO: failed to find NML: ", nml)

  end subroutine inputs_seek_to

end module bl_IO_module

