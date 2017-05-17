
subroutine f(data, lo, hi) bind(c)
  use amrex_fort_module, only : amrex_spacedim, amrex_real
  integer, intent(in) :: lo(amrex_spacedim), hi(amrex_spacedim)
  real(amrex_real) :: data(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

  print *, "Fortran: space dim ", amrex_spacedim
  print *, "Fortran: Fab lo and hi ", lo, hi
  print *, "Fortran: A slice of data ", data(lo(1),lo(2),:)
end subroutine f
