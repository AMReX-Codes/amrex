
module amrex_error_module

  use iso_c_binding
  use amrex_string_module
  
  implicit none

  interface amrex_error
     module procedure amrex_error0
     module procedure amrex_error1_ch
     module procedure amrex_error1_i
     module procedure amrex_error1_r
  end interface amrex_error

  private

  public :: amrex_error, amrex_abort, amrex_warning

  interface
     subroutine amrex_fi_error (message) bind(c)
       import
       character(kind=c_char), intent(in) :: message(*)
     end subroutine amrex_fi_error

     subroutine amrex_fi_abort (message) bind(c)
       import
       character(kind=c_char), intent(in) :: message(*)
     end subroutine amrex_fi_abort       

     subroutine amrex_fi_warning (message) bind(c)
       import
       character(kind=c_char), intent(in) :: message(*)
     end subroutine amrex_fi_warning
  end interface

contains

  subroutine amrex_error0 (message)
    character(len=*), intent(in) :: message
    call amrex_fi_error(amrex_string_f_to_c(message))
  end subroutine amrex_error0

  subroutine amrex_error1_ch (message, str)
    character(len=*), intent(in) :: message, str
    call amrex_fi_error(amrex_string_f_to_c(message // " " // trim(str)))
  end subroutine amrex_error1_ch

  subroutine amrex_error1_i (message, i)
    character(len=*), intent(in) :: message
    integer, intent(in) :: i
    character(len=16) :: imsg
    write(imsg,*) i
    call amrex_fi_error(amrex_string_f_to_c(message // " " // trim(imsg)) )
  end subroutine amrex_error1_i

  subroutine amrex_error1_r (message, r)
    use amrex_fort_module, only : amrex_real
    character(len=*), intent(in) :: message
    real(amrex_real), intent(in) :: r
    character(len=30) :: rmsg
    write(rmsg,*) r
    call amrex_fi_error(amrex_string_f_to_c(message) // " " // trim(rmsg) )
  end subroutine amrex_error1_r

  subroutine amrex_abort (message)
    character(len=*), intent(in) :: message
    call amrex_fi_abort(amrex_string_f_to_c(message))
  end subroutine amrex_abort

  subroutine amrex_warning (message)
    character(len=*), intent(in) :: message
    call amrex_fi_warning(amrex_string_f_to_c(message))
  end subroutine amrex_warning

end module amrex_error_module

! For backward compatibility
subroutine bl_error (message)
  use amrex_error_module
  implicit none
  character(len=*), intent(in) :: message
  call amrex_error(message)
end subroutine bl_error

subroutine bl_abort (message)
  use amrex_error_module
  implicit none
  character(len=*), intent(in) :: message
  call amrex_abort(message)
end subroutine bl_abort

subroutine bl_warning (message)
  use amrex_error_module
  implicit none
  character(len=*), intent(in) :: message
  call amrex_warning(message)
end subroutine bl_warning
