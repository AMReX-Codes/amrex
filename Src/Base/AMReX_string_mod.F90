
module amrex_string_module

  use iso_c_binding

  implicit none

  private

  character(c_char), public, parameter :: amrex_c_null_char_array(1) = (/ c_null_char /)

  public :: amrex_string_f_to_c, amrex_string_c_to_f

contains

  function amrex_string_f_to_c (fstr) result(cstr)
    character(*), intent(in) :: fstr
    character(c_char) :: cstr(len_trim(fstr)+1)
    integer :: i, n
    n = len_trim(fstr)
    do i = 1, n
       cstr(i) = fstr(i:i)
    end do
    cstr(n+1) = c_null_char
  end function amrex_string_f_to_c

  function amrex_string_c_to_f (cstr) result(fstr)
    character(c_char), intent(in) :: cstr(:)
    character(len=size(cstr)-1) :: fstr
    integer :: i, n
    n = size(cstr)-1   ! skip the null character
    fstr = ""
    do i = 1, n
       if (cstr(i) == c_null_char) exit
       fstr(i:i) = transfer(cstr(i), fstr)
    enddo
  end function amrex_string_c_to_f

end module amrex_string_module
