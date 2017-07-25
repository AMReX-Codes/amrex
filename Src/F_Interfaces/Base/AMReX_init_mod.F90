
module amrex_init_module

  use iso_c_binding
  implicit none

  private
  public :: amrex_init, amrex_finalize

contains

  subroutine amrex_init ()
    use amrex_string_module
    interface
       subroutine amrex_fi_init(cmd) bind(c)
         import
         implicit none
         character(c_char), intent(in) :: cmd(*)
       end subroutine amrex_fi_init
    end interface
    integer :: len_cmd, status_cmd
    character(len=1024) :: cmd
    type(amrex_string) :: cmd_string
    CALL get_command(cmd, len_cmd, status_cmd)
    if (status_cmd .ne. 0) then
       stop "amrex_init: get_command failed"
    else
       call amrex_string_build(cmd_string, cmd)
       call amrex_fi_init(cmd_string%data)
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
  end subroutine amrex_finalize

end module amrex_init_module
