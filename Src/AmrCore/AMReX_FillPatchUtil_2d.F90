
subroutine amrex_interp_div_free_bfield (lo, hi, bx, bxlo, bxhi, by, bylo, byhi, &
     cx, cxlo, cxhi, cy, cylo, cyhi, dx, rr) bind(c)
  use amrex_fort_module, only : amrex_real
  implicit none

  integer, intent(in) :: lo(2), hi(2), bxlo(2), bxhi(2), bylo(2), byhi(2), &
       cxlo(2), cxhi(2), cylo(2), cyhi(2), rr
  real(amrex_real), intent(in) :: dx(2)
  real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
  real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2))
  real(amrex_real), intent(in   ) :: cx(cxlo(1):cxhi(1),cxlo(2):cxhi(2))
  real(amrex_real), intent(in   ) :: cy(cylo(1):cyhi(1),cylo(2):cyhi(2))

  integer :: i,j

end subroutine amrex_interp_div_free_bfield
