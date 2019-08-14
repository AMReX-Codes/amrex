
subroutine my_func() bind(c)

  use amrex_fort_module, only : amrex_spacedim, amrex_real

  print *, "Fortran: space dim ", amrex_spacedim

end subroutine my_func
