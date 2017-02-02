
subroutine amrex_fmain () bind(c)

  use amrex_module
  use my_amr_module

  implicit none

  type(my_amr) :: amr

  call amr%build()

  call amr%destroy()

end subroutine amrex_fmain
