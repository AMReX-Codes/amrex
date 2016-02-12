
module backtrace_module
  
  use iso_c_binding
  use parallel

  implicit none

  private

  public :: backtrace_init, abort_fortranboxlib, fpe_trap, quiet_nan

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

  subroutine abort_fortranboxlib () bind(c,name="abort_fortranboxlib")
    call parallel_abort()
  end subroutine abort_fortranboxlib

  function fpe_trap () result(r)
    interface
       integer(c_int) function get_fpe_trap () bind(c)
         use iso_c_binding, only : c_int
       end function get_fpe_trap
    end interface
    logical :: r
    integer(c_int) flags
    flags = get_fpe_trap()
    r = flags .ne. 0
  end function fpe_trap

  function quiet_nan () result(q)
    interface
       real(c_double) function get_quiet_nan () bind(c)
         use iso_c_binding, only : c_double
       end function get_quiet_nan
    end interface
    double precision :: q
    q = get_quiet_nan()
  end function quiet_nan

end module backtrace_module

