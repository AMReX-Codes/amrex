
module amrex_error_module

  use iso_c_binding
  use amrex_string_module
  
  implicit none

  private

  public :: amrex_error, amrex_abort

  interface
     subroutine amrex_fi_error (message) bind(c)
       import
       character(kind=c_char), intent(in) :: message(*)
     end subroutine amrex_fi_error

     subroutine amrex_fi_abort (message) bind(c)
       import
       character(kind=c_char), intent(in) :: message(*)
     end subroutine amrex_fi_abort       
  end interface

contains

  subroutine amrex_error (message)
    character(len=*), intent(in) :: message
    call amrex_fi_error(amrex_string_f_to_c(message))
  end subroutine amrex_error

  subroutine amrex_abort (message)
    character(len=*), intent(in) :: message
    call amrex_fi_abort(amrex_string_f_to_c(message))
  end subroutine amrex_abort

end module amrex_error_module
