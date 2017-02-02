
subroutine amrex_fmain () bind(c)

  use amrex_module
  use amrex_famrcore_module
  use my_amr_module

  implicit none

  type(amrex_mfiter) :: mfi, mfi2

  call amrex_famrcore_init()

  call my_amr_init()

  mfi = mfi2

  call my_amr_finalize()

  call amrex_famrcore_finalize()

end subroutine amrex_fmain
