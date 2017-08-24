module cns_eb_diff_wall_module
  use amrex_fort_module, only : rt=>amrex_real
  use cns_module, only : qvar, qu, qv, qw, qtemp
  implicit none
  private
  public :: compute_diff_wallflux

contains

  subroutine compute_diff_wallflux !(divw, q, qlo, qhi, &
!       ax, axlo, axhi, ay, aylo, ayhi, az, azlo, azhi)

  end subroutine compute_diff_wallflux

end module cns_eb_diff_wall_module
