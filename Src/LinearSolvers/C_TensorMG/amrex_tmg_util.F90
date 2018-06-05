
subroutine amrex_tmg_polyInterpCoeff (xInt, x, N, c)
  use amrex_fort_module, only : amrex_real
  use amrex_lo_util_module, only : polyInterpCoeff
  implicit none
  integer :: N
  real(amrex_real) :: xInt, x(N), c(N)
  call polyInterpCoeff(xInt, x, N, c)
end subroutine amrex_tmg_polyInterpCoeff
