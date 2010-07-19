!! Module for string/character manipulations
module bl_string_module

  implicit none

  integer, parameter, private :: EOS = -1

contains

  !! Converts an integer encoding of a character string
  !! to a Fortran string
  subroutine int2str(str, iarr, n)
    use bl_error_module
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
    use bl_error_module
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
  pure function to_lower (c) result(r)
    character :: r
    character, intent(in) :: c
    r = c
    if ( is_upper(c) ) then
       r = char ( ichar(c) + 32 )
    end if
  end function to_lower

  !! Converts character to uppercase
  pure function to_upper (c) result(r)
    character :: r
    character, intent(in) :: c
    r = c
    if ( is_lower(c) ) then
       r = char ( ichar(c) - 32 )
    end if
  end function to_upper

  !! Checks for an uppercase letter.
  pure function is_upper(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = (ichar(c) >= ichar('A') .and. ichar(c) <= ichar('Z'))
  end function is_upper

  !! Checks for a lowercase letter.
  pure function is_lower(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = (ichar(c) >= ichar('a') .and. ichar(c) <= ichar('z'))
  end function is_lower

  !! Checks for a digit (0 through 9)
  pure function is_digit(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = (ichar(c) >= ichar('0') .and. ichar(c) <= ichar('9'))
  end function is_digit

  !! Checks for an alphabetic character
  pure function is_alpha(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = is_lower(c) .or. is_upper(c)
  end function is_alpha

  !! Checks for an alphanumceric character
  pure function is_alnum(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = is_alpha(c) .or. is_digit(c)
  end function is_alnum

  !! Checks for a blank character
  pure function is_blank(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = ( c == ' ' .or. ichar(c) == 9 )
  end function is_blank

  !! Checks for white-space characters.
  pure function is_space(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = ( c == ' ' .or. (ichar(c) >= 9 .and. ichar(c) < 13) )
  end function is_space

  !! Checks for any printable character including space.
  pure function is_print(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = (ichar(c) >= 32 .and. ichar(c) <= 126 )
  end function is_print

  !! Checks for any printable character except space.
  pure function is_graph(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = (ichar(c) > 32 .and. ichar(c) <= 126 )
  end function is_graph

  !! Checks fora ny printable character that is not a space
  !! or an alphanumeric character
  pure function is_punct(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = is_print(c) .and. .not. is_space(c) .and. .not. is_alnum(c)
  end function is_punct

  !! Checks for a hexadecima digit.
  pure function is_xdigit(c) result(r)
    logical :: r
    character, intent(in) :: c
    r = is_digit(c) .or. &
         ( &
         ichar('a') <= ichar(to_lower(c)) .and. &
         ichar(to_lower(c)) <= ichar('f') &
         )
  end function is_xdigit

  !! case insensitive string comparison
  pure function eq_i ( s1, s2 ) result(r)
    logical r
    character(len=*), intent(in) :: s1, s2
    integer :: i, l1, l2, lc

    l1 = len_trim (s1)
    l2 = len_trim (s2)
    r = .false.
    if ( l1 /= l2 ) return
    lc = l1
    do i = 1, lc
       if ( to_upper(s1(i:i)) /= to_upper(s2(i:i)) ) then
          return
       end if
    end do
    r = .true.

  end function eq_i

end module bl_string_module

