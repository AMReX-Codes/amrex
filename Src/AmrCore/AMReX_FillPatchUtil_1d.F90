
subroutine amrex_interp_div_free_bfield (lo, hi, bx, bxlo, bxhi, cx, cxlo, cxhi, dx, rr) bind(c)
  use amrex_fort_module, only : amrex_real
  implicit none

  integer, intent(in) :: lo(1), hi(1), bxlo(1), bxhi(1), cxlo(1), cxhi(1), rr
  real(amrex_real), intent(in) :: dx(1)
  real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1))
  real(amrex_real), intent(in   ) :: cx(cxlo(1):cxhi(1))

  integer :: i

  stop "amrex_interp_div_free_bfield in 1d not implemented"

end subroutine amrex_interp_div_free_bfield
