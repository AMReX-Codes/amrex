module bl_system_module

  implicit none

  integer, parameter :: BL_CWD_SIZE = 256

contains

  subroutine get_cwd(cwd)

    use bl_string_module
    use bl_error_module

    character (len=*), INTENT(inout) :: cwd

    integer :: icwd(BL_CWD_SIZE)

    if (len(cwd) /= BL_CWD_SIZE) then
       call bl_error("ERROR: get_cwd argument cwd wrong size")
    endif

    call get_present_dir(icwd, BL_CWD_SIZE)

    call int2str(cwd, icwd, BL_CWD_SIZE)

  end subroutine get_cwd

end module bl_system_module

