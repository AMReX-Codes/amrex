
module amrex_init_module

  use iso_c_binding
  use amrex_parallel_module
  implicit none

  private
  public :: amrex_init, amrex_finalize, amrex_initialized

  logical, save, private :: initialized = .false.

contains

  subroutine amrex_init (comm, arg_parmparse)
    use amrex_string_module
    interface
       subroutine amrex_fi_init(cmd, fcomm, argpp) bind(c)
         import
         implicit none
         character(c_char), intent(in) :: cmd(*)
         integer, intent(in), value :: fcomm, argpp
       end subroutine amrex_fi_init
    end interface
    integer, optional, intent(in) :: comm
    logical, optional, intent(in) :: arg_parmparse
    integer :: len_cmd, status_cmd, l_arg_parmparse
    character(len=1024) :: cmd
    type(amrex_string) :: cmd_string
    call amrex_parallel_init(comm)
    call get_command(cmd, len_cmd, status_cmd)
    if (status_cmd .ne. 0) then
       stop "amrex_init: get_command failed"
    else
       call amrex_string_build(cmd_string, cmd)
       l_arg_parmparse = 1
       if (present(arg_parmparse)) then
          if (arg_parmparse) then
             l_arg_parmparse = 1
          else
             l_arg_parmparse = 0
          end if
       end if
       call amrex_fi_init(cmd_string%data, amrex_parallel_communicator(), l_arg_parmparse)
       initialized = .true.
    end if
  end subroutine amrex_init
  
  subroutine amrex_finalize ()
    interface
       subroutine amrex_fi_finalize() bind(c)
         import
         implicit none
       end subroutine amrex_fi_finalize
    end interface
    call amrex_fi_finalize()
    call amrex_parallel_finalize()
    initialized = .false.
  end subroutine amrex_finalize

  pure logical function amrex_initialized ()
    amrex_initialized = initialized
  end function amrex_initialized

end module amrex_init_module
