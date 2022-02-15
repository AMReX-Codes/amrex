
program main

  use amrex_base_module, only : amrex_init, amrex_finalize

  use mytest_module, only : init, finalize, solve, write_plotfile

  implicit none

  call amrex_init()

  call init()

  call solve()

  call write_plotfile()

  call finalize()

  call amrex_finalize()

end program main
