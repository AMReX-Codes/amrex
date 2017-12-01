
module amrex_string_module

  use iso_c_binding

  implicit none

  private

  character(kind=c_char), public, parameter :: amrex_c_null_char_array(1) = (/ c_null_char /)

  type, public :: amrex_string
     character(kind=c_char), allocatable :: data(:)
  end type amrex_string

  public :: amrex_string_f_to_c, amrex_string_c_to_f, amrex_string_build, amrex_string_f_to_lower

contains

  function amrex_string_f_to_c (fstr) result(cstr)
    character(*), intent(in) :: fstr
    character(kind=c_char) :: cstr(len_trim(fstr)+1)
    integer :: i, n
    n = len_trim(fstr)
    do i = 1, n
       cstr(i) = fstr(i:i)
    end do
    cstr(n+1) = c_null_char
  end function amrex_string_f_to_c

  function amrex_string_c_to_f (cstr) result(fstr)
    character(kind=c_char), intent(in) :: cstr(:)
    character(len=size(cstr)-1) :: fstr
    integer :: i, n
    n = size(cstr)-1   ! skip the null character
    fstr = ""
    do i = 1, n
       if (cstr(i) == c_null_char) exit
       fstr(i:i) = cstr(i)
    enddo
  end function amrex_string_c_to_f

  subroutine amrex_string_build (str, fstr)
    type(amrex_string), intent(inout) :: str
    character(*), intent(in) :: fstr
    integer :: i, n
    if (allocated(str%data)) deallocate(str%data)
    n = len_trim(fstr)
    allocate(str%data(n+1))
    do i = 1, n
       str%data(i) = fstr(i:i)
    end do
    str%data(n+1) = c_null_char
  end subroutine amrex_string_build


  function amrex_string_f_to_lower (fstr) result(flo)
    character(len=*), intent(in) :: fstr
    character(len=len_trim(fstr)) :: flo
    integer :: i, n
    n = len_trim(fstr)
    do i = 1, n
       select case(fstr(i:i))
          case ("A":"Z")
             flo(i:i) = achar(iachar(fstr(i:i))+32)
          case default
             flo(i:i) = fstr(i:i)
       end select
    end do
  end function amrex_string_f_to_lower

end module amrex_string_module
