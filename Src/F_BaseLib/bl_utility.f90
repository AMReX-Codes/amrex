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

!! Module for string/character manipulations
module bl_string_module

  use bl_error_module

  implicit none

  integer, parameter, private :: EOS = -1

contains

  !! Converts an integer encoding of a character string
  !! to a Fortran string
  subroutine int2str(str, iarr, n)
    character(len=*), intent(out) :: str
    integer, intent(in) :: n
    integer, intent(in) :: iarr(n)
    integer :: i

    if ( len(str) < n ) then
       call bl_error("INT2STR: iarr to large for str: len = ", len(str))
    end if
    do i = 1, n
       if ( iarr(i) == EOS ) exit
       str(i:i) = char(iarr(i))
    end do

  end subroutine int2str

  !! Converts a Fortran string to an integer encoding.
  subroutine str2int(iarr, n, str)
    character(len=*), intent(in) :: str
    integer, intent(in) :: n
    integer :: i, j
    integer, intent(out) :: iarr(n)

    if ( n <= len_trim(str) ) then
       call bl_error("STR2INT: str to large for iarr: size(iarr) = ", n)
    end if

    iarr = 0
    j = 1
    do i = 1, len_trim(str)
       iarr(j) = ichar(str(i:i))
       j = j + 1
    end do

    iarr(j) = EOS

  end subroutine str2int

  !! Converts character to lowercase
  function to_lower ( c ) result(r)
    character :: r
    character, intent(in) :: c
    r = c
    if ( is_upper(c) ) then
       r = char ( ichar(c) + 32 )
    end if
  end function to_lower

  !! Converts character to uppercase
  function to_upper ( c ) result(r)
    character r
    character, intent(in) :: c
    r = c
    if ( is_lower(c) ) then
       r = char ( ichar(c) - 32 )
    end if
  end function to_upper

  !! Checks for an uppercase letter.
  function is_upper(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = (ichar(c) >= ichar('A') .and. ichar(c) <= ichar('Z'))
  end function is_upper

  !! Checks for a lowercase letter.
  function is_lower(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = (ichar(c) >= ichar('a') .and. ichar(c) <= ichar('z'))
  end function is_lower

  !! Checks for a digit (0 through 9)
  function is_digit(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = (ichar(c) >= ichar('0') .and. ichar(c) <= ichar('9'))
  end function is_digit

  !! Checks for an alphabetic character
  function is_alpha(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = is_lower(c) .or. is_upper(c)
  end function is_alpha

  !! Checks for an alphanumceric character
  function is_alnum(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = is_alpha(c) .or. is_digit(c)
  end function is_alnum

  !! Checks for a blank character
  function is_blank(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = ( c == ' ' .or. ichar(c) == 9 )
  end function is_blank

  !! Checks for white-space characters.
  function is_space(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = ( c == ' ' .or. (ichar(c) >= 9 .and. ichar(c) < 13) )
  end function is_space

  !! Checks for any printable character including space.
  function is_print(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = (ichar(c) >= 32 .and. ichar(c) <= 126 )
  end function is_print

  !! Checks for any printable character except space.
  function is_graph(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = (ichar(c) > 32 .and. ichar(c) <= 126 )
  end function is_graph

  !! Checks fora ny printable character that is not a space
  !! or an alphanumeric character
  function is_punct(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = is_print(c) .and. .not. is_space(c) .and. .not. is_alnum(c)
  end function is_punct

  !! Checks for a hexadecima digit.
  function is_xdigit(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = is_digit(c) .or. &
         ( &
         ichar('a') <= ichar(to_lower(c)) .and. &
         ichar(to_lower(c)) <= ichar('f') &
         )
  end function is_xdigit

  !! case insensitive string comparison
  function eq_i ( s1, s2 ) result(r)
    logical r
    character(len=*), intent(in) :: s1, s2
    integer i, l1, l2, lc

    l1 = len ( s1 )
    l2 = len ( s2 )
    lc = min ( l1, l2 )
    r = .false.
    do i = 1, lc
       if ( to_upper(s1(i:i)) /= to_upper(s2(i:i)) ) then
          return
       end if
    end do
    do i = lc + 1, l1
       if ( s1(i:i) /= ' ' ) then
          return
       end if
    end do
    do i = lc + 1, l2
       if ( s2(i:i) /= ' ' ) then
          return
       end if
    end do
    r = .true.

  end function eq_i

end module bl_string_module

!! Presents a model of stream I/O
!!
!! A unit consists of whitespace separated characters that can
!! be scanned.  Single characters can be put back, enabling
!! primitive parsing.

module bl_stream_module
  
  use bl_error_module
  use bl_string_module
  use bl_IO_module

  implicit none

  type bl_stream
     integer   :: unit = -1
     character ::  pb  = achar(0)
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

  function bl_stream_the_unit(strm) result(r)
    integer :: r
    type(bl_stream), intent(in) :: strm
    r = strm%unit
  end function bl_stream_the_unit

  subroutine bl_stream_eat_whitespace(strm)
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
100 end do
  end subroutine bl_stream_eat_whitespace

  function bl_stream_peek_chr(strm, set) result(r)
    character :: r
    type(bl_stream), intent(inout) :: strm
    character(len=*), intent(in), optional :: set
    r = bl_stream_scan_chr(strm, set)
    call bl_stream_putback_chr(strm, r)
  end function bl_stream_peek_chr

  function bl_stream_scan_int(strm) result(r)
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
          if ( .not. started ) then
             call bl_error("BL_STREAM_SCAN_INT: first character in scan: ", c)
          end if
          call bl_stream_putback_chr(strm, c)
          exit
       end if
       started = .true.
       r = r*10 + (ichar(c)-ichar('0'))
    end do
100 continue
  end function bl_stream_scan_int

  function bl_stream_scan_chr(strm, set) result(r)
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
100 end do
  end function bl_stream_scan_chr

  subroutine bl_stream_expect_str(strm, str)
    type(bl_stream), intent(inout) :: strm
    character(len=*), intent(in) :: str
    character :: c
    integer :: i
    c = bl_stream_scan_chr(strm)
    do i = 1, len(str)
       if ( c /= str(i:i) ) then
          call bl_error('expect_chr: expected "'// str// '" got: ', &
               '"'//c//'"')
        end if
       if ( i == len(str) ) exit
       c = bl_stream_scan_chr(strm)
    end do
  end subroutine bl_stream_expect_str

  subroutine bl_stream_expect_chr_v(strm, chars)
    type(bl_stream), intent(inout) :: strm
    character, intent(in) :: chars(:)
    integer i
    do i = 1, len(chars)
       call bl_stream_expect_str(strm, chars(i))
    end do
  end subroutine bl_stream_expect_chr_v

  subroutine bl_stream_putback_chr(strm, char)
    type(bl_stream), intent(inout) :: strm
    character, intent(in) :: char

    if ( strm%lpb ) then
       call bl_error('PUTBACK_CHR: already pushed back')
    end if
    strm%pb = char
    strm%lpb = .TRUE.

  end subroutine bl_stream_putback_chr

end module bl_stream_module

module bl_mem_stat_module
  
  use bl_types
  use parallel

  implicit none

  type mem_stats
     integer(kind =ll_t) :: num_alloc = 0_ll_t
     integer(kind =ll_t) :: num_dealloc = 0_ll_t
     integer(kind =ll_t) :: cnt_alloc = 0_ll_t
     integer(kind =ll_t) :: cnt_dealloc = 0_ll_t
  end type mem_stats

  interface print
     module procedure mem_stats_print
  end interface

  interface mem_stats_alloc
     module procedure mem_stats_alloc_ll
     module procedure mem_stats_alloc_i
     module procedure mem_stats_alloc_c
  end interface

  interface mem_stats_dealloc
     module procedure mem_stats_dealloc_ll
     module procedure mem_stats_dealloc_i
     module procedure mem_stats_dealloc_c
  end interface

contains
!!345678901234567890123456789012345678901234567890123456789012345678901234567890
!!       1         2         3         4         5         6         7         8

  subroutine mem_stats_print(ms, str, unit, advance, total)
    use bl_IO_module
    use bl_string_module
    type(mem_stats), intent(in) :: ms
    character(len=*), intent(in), optional :: str
    character(len=*), intent(in), optional :: advance
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: total
    integer :: un
    character(len=3) adv
    logical :: ltotal
    ltotal = .true.; if ( present(total) ) ltotal = total
    if ( parallel_IOProcessor() ) then
       un = unit_stdout(unit)
       adv = unit_advance(advance)
       if ( present(str) ) then
          write(unit = un, fmt = '(A,": ")', advance = 'no') str
       end if
       if ( ltotal ) then
          write(unit=un, fmt='(4i15)', advance = 'no') &
               ms%num_alloc,                           &
               ms%num_alloc-ms%num_dealloc,            &
               ms%cnt_alloc,                           &
               ms%cnt_alloc-ms%cnt_dealloc
       else
          write(unit=un, fmt='(2i15)', advance = 'no') &
               ms%num_alloc-ms%num_dealloc,            &
               ms%cnt_alloc-ms%cnt_dealloc
       end if
       if ( eq_i(adv, "YES") ) write(unit=un, fmt='()')
    end if
  end subroutine mem_stats_print

  subroutine mem_stats_alloc_ll(ms, vol)
    type(mem_stats), intent(inout) :: ms
    integer(kind=ll_t) :: vol
    ms%num_alloc = ms%num_alloc + vol
    ms%cnt_alloc = ms%cnt_alloc + 1
  end subroutine mem_stats_alloc_ll

  subroutine mem_stats_alloc_i(ms, vol)
    type(mem_stats), intent(inout) :: ms
    integer :: vol
    ms%num_alloc = ms%num_alloc + int(vol,kind=ll_t)
    ms%cnt_alloc = ms%cnt_alloc + 1
  end subroutine mem_stats_alloc_i

  subroutine mem_stats_alloc_c(ms)
    type(mem_stats), intent(inout) :: ms
    ms%cnt_alloc = ms%cnt_alloc + 1
  end subroutine mem_stats_alloc_c

  subroutine mem_stats_dealloc_ll(ms, vol)
    type(mem_stats), intent(inout) :: ms
    integer(kind=ll_t) :: vol
    ms%num_dealloc = ms%num_dealloc + vol
    ms%cnt_dealloc = ms%cnt_dealloc + 1
  end subroutine mem_stats_dealloc_ll

  subroutine mem_stats_dealloc_i(ms, vol)
    type(mem_stats), intent(inout) :: ms
    integer :: vol
    ms%num_dealloc = ms%num_dealloc + int(vol,kind=ll_t)
    ms%cnt_dealloc = ms%cnt_dealloc + 1
  end subroutine mem_stats_dealloc_i

  subroutine mem_stats_dealloc_c(ms)
    type(mem_stats), intent(inout) :: ms
    ms%cnt_dealloc = ms%cnt_dealloc + 1
  end subroutine mem_stats_dealloc_c

end module bl_mem_stat_module

!! A module for timers
module bl_timer_module

  use bl_types
  use bl_error_module
  use parallel

  implicit none

  integer, parameter :: TIMER_NAME_MAX_LEN = 20
  real(kind=dp_t), private, parameter :: zero = 0.0_dp_t
  real(kind=dp_t), private, parameter :: MIL_SEC = 1.0e3_dp_t

  type timer
     character(len=TIMER_NAME_MAX_LEN) :: name = ""
     real(kind=dp_t) :: strt = zero
     real(kind=dp_t) :: time = zero
     real(kind=dp_t) :: cum_time = zero
     integer(kind=ll_t) :: cnt = 0_ll_t
     logical :: running = .FALSE.
  end type timer

  interface
     subroutine cpu_second(s)
       use bl_types
       real(kind=dp_t) :: s
     end subroutine cpu_second
     subroutine wall_second(s)
       use bl_types
       real(kind=dp_t) :: s
     end subroutine wall_second
     subroutine cpu_second_tick(s)
       use bl_types
       real(kind=dp_t) :: s
     end subroutine cpu_second_tick
     subroutine wall_second_tick(s)
       use bl_types
       real(kind=dp_t) :: s
     end subroutine wall_second_tick
  end interface

  private :: cpu_second, wall_second
 
contains

  function timer_construct(str) result(r)
    character(len=*), intent(in) :: str
    type(timer) :: r
    r%name = str
    r%time = zero
    r%strt = zero
    r%cum_time = zero
    r%cnt = 0_ll_t
    r%running = .FALSE.
  end function timer_construct

  subroutine timer_start(tm)
    type(timer), intent(inout) :: tm

    if ( tm%running ) &
         call bl_error("TIMER_START: should not be be running before start")
    call wall_second(tm%strt)
    tm%running = .TRUE.
    tm%time = zero
  end subroutine timer_start

  subroutine timer_stop(tm)
    type(timer), intent(inout) :: tm

    if ( .not. tm%running ) &
         call bl_error("TIMER_STOP: should be running before stop")
    call wall_second(tm%time)
    tm%time = tm%time-tm%strt
    tm%strt = zero
    tm%running = .FALSE.
    tm%cum_time = tm%cum_time + tm%time
    tm%cnt = tm%cnt + 1_ll_t

  end subroutine timer_stop

  ! This is dangerous since it forces a reduction!
  function timer_value(tm, total) result(r)
    real(kind=dp_t) :: r
    type(timer), intent(in) :: tm
    logical, intent(in), optional :: total
    logical ltotal
    real(kind=dp_t) :: r1
    ltotal = .FALSE.; if ( present(total) ) ltotal = total
    if ( tm%running ) &
         call bl_error("TIMER_VALUE: should be stopped before getting value")
    if ( ltotal ) then
       r1 = tm%cum_time
    else
       r1 = tm%time
    end if
    call parallel_reduce(r, r1, MPI_MAX)
    r = r1*MIL_SEC
  end function timer_value

  subroutine timer_print(tm, str, unit, advance, total)
    use bl_IO_module
    use bl_string_module
    type(timer), intent(in) :: tm
    character(len=*), intent(in), optional :: str
    character(len=*), intent(in), optional :: advance
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: total
    integer un
    character(len=3) adv
    logical :: ltotal
    real :: r

    !! Timer allows 9.9999999e9 as largest number of milliseconds print
    !! this is 3 years.
    !!      3301.5550
    !! F15.4 ( I have used a few timers with resolutions as low as .95 us.
    !!> 1234567890
    !!            1
    !!<            2345
    
    r = timer_value(tm, total)
    if ( parallel_IOProcessor() ) then
       un = unit_stdout(unit)
       adv = unit_advance(advance)
       ltotal = .false.; if ( present(total) ) ltotal = total
       if ( present(str) ) then
          write(unit=un, fmt='("(*", A, "*)")', advance = adv) str
       end if
       if ( ltotal ) then
          write(unit=un, &
               fmt='("timer[",A,",",i15,",",F15.4,"]")', advance = 'no') &
               adjustr(tm%name), tm%cnt, tm%cum_time*MIL_SEC
       else
          write(unit=un, fmt='("timer[",A,",",F15.4,"]")', advance = 'no') &
               adjustr(tm%name), tm%time*MIL_SEC
       end if
       if ( eq_i(adv, 'YES') ) write(unit = un, fmt='()')
    end if

  end subroutine timer_print

  function timer_tick() result(r)
    real(kind=dp_t) :: r
    call wall_second_tick(r)
    r = r*MIL_SEC
  end function timer_tick

end module bl_timer_module

module kiss_module

  ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
  ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
  ! (2) A 3-shift shift-register generator, period 2^32-1,
  ! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497 > 2^59
  !  Overall period > 2^123;  Default seeds x,y,z,w.
  !  Set your own seeds with statement i=kisset(ix,iy,iz,iw).

  implicit none

  private
  public :: kiss, kiss_init, kiss_rn

  type :: kiss_rn
     private
     integer :: x = 123456789
     integer :: y = 362436069
     integer :: z = 521288629
     integer :: w = 916191069
  end type kiss_rn
  type(kiss_rn), save :: the_kiss

  interface kiss
     module procedure kiss_d
     module procedure kiss_k
  end interface

  interface kiss_init
     module procedure kiss_init_d
     module procedure kiss_init_k
  end interface

  interface kiss_fill
     module procedure kiss_fill_0d
     module procedure kiss_fill_1d
     module procedure kiss_fill_2d
     module procedure kiss_fill_3d
     module procedure kiss_fill_4d
     module procedure kiss_fill_5d
     module procedure kiss_fill_6d
     module procedure kiss_fill_7d
  end interface

contains

  function kiss_k (k, range) result(r)
    type(kiss_rn), intent(inout) :: k
    integer, intent(in), optional :: range
    integer :: r
    k%x = 69069 * k%x + 1327217885
    k%y = m (m (m (k%y, 13), - 17), 5)
    k%z = 18000 * iand (k%z, 65535) + ishft (k%z, - 16)
    k%w = 30903 * iand (k%w, 65535) + ishft (k%w, - 16)
    r = k%x + k%y + ishft (k%z, 16) + k%w
    if ( present(range) ) r = mod(r,range)
  contains
    function m(k, n)
      integer :: m
      integer, intent(in) :: k, n
      m = ieor (k, ishft (k, n) )
    end function m
  END FUNCTION kiss_k

  function kiss_d (range) result(r)
    integer :: r
    integer, intent(in), optional :: range
    r = kiss_k(the_kiss, range)
  end function kiss_d

  subroutine kiss_init_k (k, ix, iy, iz, iw)
    type(kiss_rn), intent(out) :: k
    integer, intent(in) :: ix, iy, iz, iw
    k%x = ix
    k%y = iy
    k%z = iz
    k%w = iw
  end subroutine kiss_init_k

  subroutine kiss_init_d (ix, iy, iz, iw)
    integer, intent(in) :: ix, iy, iz, iw
    call kiss_init_k(the_kiss, ix, iy, iz, iw)
  end subroutine kiss_init_d

  subroutine kiss_fill_0d(a, range)
    integer, intent(out) :: a
    integer, intent(in), optional :: range
    a = kiss(range)
  end subroutine kiss_fill_0d

  subroutine kiss_fill_1d(a, range)
    integer, intent(out)  :: a(:)
    integer, intent(in), optional :: range
    integer :: i
    do i = 1, size(a,1)
       a(i) = kiss(range)
    end do

  end subroutine kiss_fill_1d

  subroutine kiss_fill_2d(a, range)
    integer, intent(out)  :: a(:,:)
    integer, intent(in), optional :: range
    integer :: i, j
    do j = 1, size(a,2)
       do i = 1, size(a,1)
          a(i,j) = kiss(range)
       end do
    end do

  end subroutine kiss_fill_2d

  subroutine kiss_fill_3d(a, range)
    integer, intent(out)  :: a(:,:,:)
    integer, intent(in), optional :: range
    integer :: i, j, k
    do k = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             a(i,j,k) = kiss(range)
          end do
       end do
    end do

  end subroutine kiss_fill_3d

  subroutine kiss_fill_4d(a, range)
    integer, intent(out)  :: a(:,:,:,:)
    integer, intent(in), optional :: range
    integer :: i, j, k, l

    do l = 1, size(a,4)
       do k = 1, size(a,3)
          do j = 1, size(a,2)
             do i = 1, size(a,1)
                a(i,j,k,l) = kiss(range)
             end do
          end do
       end do
    end do

  end subroutine kiss_fill_4d

  subroutine kiss_fill_5d(a, range)
    integer, intent(out)  :: a(:,:,:,:,:)
    integer, intent(in), optional :: range
    integer :: i, j, k, l, m

    do m = 1, size(a,5)
       do l = 1, size(a,4)
          do k = 1, size(a,3)
             do j = 1, size(a,2)
                do i = 1, size(a,1)
                   a(i,j,k,l,m) = kiss(range)
                end do
             end do
          end do
       end do
    end do

  end subroutine kiss_fill_5d

  subroutine kiss_fill_6d(a, range)
    integer, intent(out)  :: a(:,:,:,:,:,:)
    integer, intent(in), optional :: range
    integer :: i, j, k, l, m, n

    do n = 1, size(a,6)
       do m = 1, size(a,5)
          do l = 1, size(a,4)
             do k = 1, size(a,3)
                do j = 1, size(a,2)
                   do i = 1, size(a,1)
                      a(i,j,k,l,m,n) = kiss(range)
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine kiss_fill_6d

  subroutine kiss_fill_7d(a, range)
    integer, intent(out)  :: a(:,:,:,:,:,:,:)
    integer, intent(in), optional :: range
    integer :: i, j, k, l, m, n, o

    do o = 1, size(a,7)
       do n = 1, size(a,6)
          do m = 1, size(a,5)
             do l = 1, size(a,4)
                do k = 1, size(a,3)
                   do j = 1, size(a,2)
                      do i = 1, size(a,1)
                         a(i,j,k,l,m,n,o) = kiss(range)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine kiss_fill_7d

end module kiss_module

module parmparse_module

  use bl_error_module

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

  function pp_arg_count(str, str_long) result(r)
    use f2kcli
    character(len=*), intent(in) :: str
    character(len=*), intent(in), optional :: str_long
    integer :: r
    integer :: narg, f
    character(len=128) :: fname
    narg = command_argument_count()
    r = 0
    do f = 1, narg
       call get_command_argument(f, value = fname)
       if ( str == fname ) then
          r = r + 1
       else if ( present(str_long) ) then
          if ( str_long == fname ) then
             r = r + 1
          end if
       end if
    end do
  end function pp_arg_count

end module parmparse_module
