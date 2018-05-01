program main

  use amrex_base_module, only : amrex_init, amrex_finalize

  use solve_helmholtz, only : solve

  use init_helmholtz, only : init, finalize, write_plotfile

  implicit none

  integer :: counter	

  call amrex_init()

  call init()

  call solve()

  call write_plotfile()

  call finalize()

  call amrex_finalize()

end program main
