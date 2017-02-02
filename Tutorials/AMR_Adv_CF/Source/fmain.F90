
subroutine amrex_fmain () bind(c)

  use amrex_module
  use amrex_famrcore_module
  use my_amr_module
  use initdata_module

  implicit none

  call amrex_famrcore_init()

  call my_amr_init()

  call initdata()

  call my_amr_finalize()

  call amrex_famrcore_finalize()

end subroutine amrex_fmain
