!! Presents a model of stream I/O
!!
!! A unit consists of whitespace separated characters that can
!! be scanned.  Single characters can be put back, enabling
!! primitive parsing.

module bl_stream_module

  implicit none

  type bl_stream
     private
     integer   :: unit = -1
     character ::  pb
     logical   :: lpb  = .False.
  end type

  character, private, parameter :: EOF = char(0)

  interface bl_stream_expect
     module procedure bl_stream_expect_str
     module procedure bl_stream_expect_chr_v
  end interface

  interface build
     module procedure bl_stream_build
  end interface

  interface destroy
     module procedure bl_stream_destroy
  end interface

contains

  subroutine bl_stream_build(strm, unit)
    use bl_IO_module
    type(bl_stream), intent(out) :: strm
    integer, intent(in), optional :: unit
    if ( present(unit) ) then
       strm%unit = unit
    else
       strm%unit = unit_new()
    end if
  end subroutine bl_stream_build

  subroutine bl_stream_destroy(strm)
    type(bl_stream), intent(inout) :: strm
    strm%unit = -1
    strm%pb   = achar(0)
    strm%lpb  = .False.
  end subroutine bl_stream_destroy

  pure function bl_stream_the_unit(strm) result(r)
    integer :: r
    type(bl_stream), intent(in) :: strm
    r = strm%unit
  end function bl_stream_the_unit

  subroutine bl_stream_eat_whitespace(strm)
    use bl_string_module
    type(bl_stream), intent(inout) :: strm
    character :: c
    do
       if ( strm%lpb ) then
          c = strm%pb
          strm%lpb = .false.
       else
          read(unit=strm%unit, fmt='(a)', advance = 'no', eor = 100) c
       end if
       if ( is_space(c) ) cycle
       call bl_stream_putback_chr(strm, c)
       exit
100    continue
    end do
  end subroutine bl_stream_eat_whitespace

  function bl_stream_peek_chr(strm, set) result(r)
    character :: r
    type(bl_stream), intent(inout) :: strm
    character(len=*), intent(in), optional :: set
    r = bl_stream_scan_chr(strm, set)
    call bl_stream_putback_chr(strm, r)
  end function bl_stream_peek_chr

  function bl_stream_scan_int(strm) result(r)
    use bl_error_module
    use bl_string_module
    integer :: r
    type(bl_stream), intent(inout) :: strm
    character :: c
    logical :: started
    call bl_stream_eat_whitespace(strm)
    r = 0
    started = .false.
    do 
       if ( strm%lpb ) then
          c = strm%pb
          strm%lpb = .false.
       else
          read(unit=strm%unit, fmt='(a)', advance = 'no', eor = 100) c
       end if
       if ( .not. is_digit(c) ) then
          if ( .not. started ) call bl_error("BL_STREAM_SCAN_INT: first character in scan: ", c)
          call bl_stream_putback_chr(strm, c)
          exit
       end if
       started = .true.
       r = r*10 + (ichar(c)-ichar('0'))
    end do
100 continue
  end function bl_stream_scan_int

  function bl_stream_scan_chr(strm, set) result(r)
    use bl_string_module
    type(bl_stream), intent(inout) :: strm
    character(len=*), intent(in), optional :: set
    character :: r
    do
       if ( strm%lpb ) then
          r = strm%pb
          strm%lpb = .FALSE.
       else
          read(unit=strm%unit, fmt='(a)', advance = 'NO', eor = 100, end = 200) r
       end if
       if ( present(set) ) then
          if ( scan(r,set) == 0 ) exit
       else
          ! skip whitespace
          if ( is_space(r) ) cycle
          exit
       end if
200    r = EOF
       exit
100    continue
    end do
  end function bl_stream_scan_chr

  subroutine bl_stream_expect_str(strm, str)
    use bl_error_module
    type(bl_stream), intent(inout) :: strm
    character(len=*), intent(in) :: str
    character :: c
    integer :: i
    c = bl_stream_scan_chr(strm)
    do i = 1, len(str)
       if ( c /= str(i:i) ) then
          call bl_error('expect_chr: expected "'// str// '" got: ', '"'//c//'"')
        end if
       if ( i == len(str) ) exit
       c = bl_stream_scan_chr(strm)
    end do
  end subroutine bl_stream_expect_str

  subroutine bl_stream_expect_chr_v(strm, chars)
    type(bl_stream), intent(inout) :: strm
    character, intent(in) :: chars(:)
    integer :: i
    do i = 1, len(chars)
       call bl_stream_expect_str(strm, chars(i))
    end do
  end subroutine bl_stream_expect_chr_v

  subroutine bl_stream_putback_chr(strm, char)
    use bl_error_module
    type(bl_stream), intent(inout) :: strm
    character, intent(in) :: char
    if ( strm%lpb ) call bl_error('PUTBACK_CHR: already pushed back')
    strm%pb = char
    strm%lpb = .TRUE.
  end subroutine bl_stream_putback_chr

end module bl_stream_module
