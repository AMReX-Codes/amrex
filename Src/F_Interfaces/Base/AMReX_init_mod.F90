
module amrex_init_module

  use iso_c_binding
  use amrex_parallel_module
  implicit none

  private
  public :: amrex_init, amrex_finalize, amrex_initialized, amrex_null_proc

  interface
     subroutine amrex_null_proc () bind(c)
     end subroutine amrex_null_proc
  end interface

  logical, save, private :: initialized = .false.

contains

  subroutine amrex_init (comm, arg_parmparse, proc_parmparse)
    use amrex_string_module
    interface
       subroutine amrex_fi_init(cmd, fcomm, argpp, procpp) bind(c)
         import
         implicit none
         character(kind=c_char), intent(in) :: cmd(*)
         integer, intent(in), value :: fcomm, argpp
         type(c_funptr), value, intent(in) :: procpp
       end subroutine amrex_fi_init
    end interface
    integer, optional, intent(in) :: comm
    logical, optional, intent(in) :: arg_parmparse
    procedure(amrex_null_proc), optional :: proc_parmparse
    integer :: len_cmd, status_cmd, l_arg_parmparse
    character(len=1024) :: cmd
    type(amrex_string) :: cmd_string
    type(c_funptr) :: cfp

    call amrex_parallel_init(comm)

    call get_command(cmd, len_cmd, status_cmd)

    call amrex_string_build(cmd_string, cmd)

    l_arg_parmparse = 1
    if (present(arg_parmparse)) then
       if (arg_parmparse) then
          l_arg_parmparse = 1
       else
          l_arg_parmparse = 0
       end if
    end if

    if (present(proc_parmparse)) then
       cfp = c_funloc(proc_parmparse)
    else
       cfp = c_null_funptr
    end if

    call amrex_fi_init(cmd_string%data, amrex_parallel_communicator(), l_arg_parmparse, cfp)
    initialized = .true.
  end subroutine amrex_init
  
  subroutine amrex_finalize ()
    use amrex_geometry_module, only: amrex_geometry_finalize
    interface
       subroutine amrex_fi_finalize() bind(c)
         import
         implicit none
       end subroutine amrex_fi_finalize
    end interface
    call amrex_fi_finalize()
    call amrex_geometry_finalize()
    call amrex_parallel_finalize()
    initialized = .false.
  end subroutine amrex_finalize

  pure logical function amrex_initialized ()
    amrex_initialized = initialized
  end function amrex_initialized

end module amrex_init_module
