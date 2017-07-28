
program main

  use amrex_base_module

  implicit none

  call amrex_init()

  if (amrex_parallel_ioprocessor()) then
     print *, "Hello world!"
  end if

  call amrex_finalize()

end program main
