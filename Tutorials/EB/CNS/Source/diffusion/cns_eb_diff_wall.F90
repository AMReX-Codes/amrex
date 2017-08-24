module cns_eb_diff_wall_module
  use amrex_fort_module, only : rt=>amrex_real
  use cns_module, only : qvar, qu, qv, qw, qtemp
  implicit none
  private
  public :: compute_diff_wallflux

contains

  subroutine compute_diff_wallflux (divw, dx, i,j,k, q, qlo, qhi, &
       lam, mu, xi, clo, chi, &
       apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi)
    integer, intent(in) :: i,j,k,qlo(3),qhi(3),clo(3),chi(3),axlo(3),axhi(3), &
         aylo(3),ayhi(3),azlo(3),azhi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(out) :: divw(5)
    real(rt), intent(in) :: q  (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qvar)
    real(rt), intent(in) :: lam(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in) :: mu (clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in) :: xi (clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(rt), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(rt), intent(in) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))

    divw = 0.d0
  end subroutine compute_diff_wallflux

end module cns_eb_diff_wall_module
