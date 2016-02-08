
module backtrace_module
  
  use iso_c_binding
  use parallel

  implicit none

  private

  public :: backtrace_init

contains

  subroutine backtrace_init ()
    interface
       subroutine set_signal_handler (ename, rank) bind(c)
         use iso_c_binding, only : c_char, c_int
         character(kind=c_char), intent(in) :: ename(*)
         integer(kind=c_int), intent(in), value :: rank
       end subroutine set_signal_handler
    end interface
    character(len=255) :: exename
    character(kind=c_char) :: exename_c(256)
    integer :: i, n
    call get_command_argument(0, exename, n)
    do i = 1, n
       exename_c(i) = exename(i:i)
    end do
    exename_c(n+1) = c_null_char
    call set_signal_handler (exename_c, parallel_myproc())
  end subroutine backtrace_init

end module backtrace_module

subroutine abort_fortranboxlib () bind(c)
  use parallel
  call parallel_abort()
end subroutine abort_fortranboxlib

