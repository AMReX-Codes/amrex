
subroutine amrex_fmain () bind(c)

  use amrex_base_module

  implicit none

  if (amrex_parallel_ioprocessor()) then
     print *, "Hello world!"
  end if

end subroutine amrex_fmain
