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

