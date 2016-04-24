
module string_module

  use iso_c_binding

  implicit none

  private

  public :: string_f_to_c

contains

  function string_f_to_c (fstr) result(cstr)
    character(*), intent(in) :: fstr
    character(c_char) :: cstr(len_trim(fstr)+1)
    integer :: i, n
    n = len_trim(fstr)
    do i = 1, n
       cstr(i) = fstr(i:i)
    end do
    cstr(n+1) = c_null_char
  end function string_f_to_c

end module string_module
