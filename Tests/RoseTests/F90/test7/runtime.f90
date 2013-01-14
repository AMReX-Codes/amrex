module runtime_init_module

  use probin_module, only : a, b, c

  implicit none

  namelist /probin/ a
  namelist /probin/ b
  namelist /probin/ c

contains

  subroutine runtime_init
    ! I/O using namelist probin
  end subroutine runtime_init

end module runtime_init_module
