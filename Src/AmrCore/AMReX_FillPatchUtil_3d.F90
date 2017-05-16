
subroutine amrex_interp_div_free_bfield (lo, hi, bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, &
     cx, cxlo, cxhi, cy, cylo, cyhi, cz, czlo, czhi) bind(c)
  use amrex_fort_module, only : amrex_real
  implicit none

  integer, intent(in) :: lo(3), hi(3), bxlo(3), bxhi(3), bylo(3), byhi(3), bzlo(3), bzhi(3), &
       cxlo(3), cxhi(3), cylo(3), cyhi(3), czlo(3), czhi(3)
  real(amrex_real), intent(inout) :: bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
  real(amrex_real), intent(inout) :: by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
  real(amrex_real), intent(inout) :: bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
  real(amrex_real), intent(in   ) :: cx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3))
  real(amrex_real), intent(in   ) :: cy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3))
  real(amrex_real), intent(in   ) :: cz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3))

  integer :: i,j,k

end subroutine amrex_interp_div_free_bfield
