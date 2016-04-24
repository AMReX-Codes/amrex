
module sting_module

  implicit none

  public

contains

  function string_f_to_c (fs) result cs
    use iso_c_binding
    character(*), intent(in) :: fs
    character(c_char) :: cs(len_trim(fs)+1)
    integer :: i, n
    n = len_trim(fs)
    do i = 1, n
       cs(i) = fs(i:i)
    end do
    cs(n+1) = c_null_char
  end function string_f_to_c

end module sting_module
