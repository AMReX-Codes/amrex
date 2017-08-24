module cns_eb_diff_wall_module
  use amrex_fort_module, only : rt=>amrex_real
  use amrex_error_module, only : amrex_abort
  use cns_module, only : qvar, qu, qv, qw, qtemp
  implicit none
  private
  public :: compute_diff_wallflux

contains

  !
  ! The wall is assumed to be adiabatic for now. Thus dTdn=0.
  ! Later we can support const T wall
  !
  ! We use no-slip boundary for velocities.
  !
  subroutine compute_diff_wallflux (divw, dx, i,j,k, q, qlo, qhi, &
       lam, mu, xi, clo, chi, bcent, blo, bhi, &
       apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi)
    integer, intent(in) :: i,j,k,qlo(3),qhi(3),clo(3),chi(3),axlo(3),axhi(3), &
         aylo(3),ayhi(3),azlo(3),azhi(3),blo(3),bhi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(out) :: divw(5)
    real(rt), intent(in) :: q  (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qvar)
    real(rt), intent(in) :: lam(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in) :: mu (clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in) :: xi (clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in) :: bcent( blo(1): bhi(1), blo(2): bhi(2), blo(3): bhi(3))
    real(rt), intent(in) :: apx  (axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(rt), intent(in) :: apy  (aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(rt), intent(in) :: apz  (azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))

    real(rt) :: apnorm, apnorminv, anrmx, anrmy, anrmz
    integer :: s

    divw = 0.d0

    apnorm = sqrt((apx(i-1,j,k)-apx(i,j,k))**2 &
         +        (apy(i,j-1,k)-apy(i,j,k))**2 &
         +        (apz(i,j,k-1)-apz(i,j,k))**2)

    if (apnorm .eq. 0.d0) then
       call amrex_abort("compute_diff_wallflux: we are in trouble.")
    end if

    apnorminv = 1.d0/apnorm
    anrmx = (apx(i-1,j,k)-apx(i,j,k)) * apnorminv  ! unit vector pointing to the wall
    anrmy = (apy(i,j-1,k)-apy(i,j,k)) * apnorminv
    anrmz = (apz(i,j,k-1)-apz(i,j,k)) * apnorminv

    ! The center of the wall
    

    if (abs(anrmx).ge.abs(anrmy) .and. abs(anrmx).ge.abs(anrmz)) then
       ! y-z plane
       s = int(sign(1.0_rt,-anrmx))
       ! y-z plane at i+s
       
       ! y-z plane at i+2s
    else if (abs(anrmy).ge.abs(anrmx) .and. abs(anrmy).ge.abs(anrmz)) then
       ! z-x plane
!       print *, i,j,k,anrmx, anrmy, anrmz
!       flush(6)
!       call amrex_abort("z-x plane")
    else
       ! x-y plane
!       print *, i,j,k,anrmx, anrmy, anrmz
!       flush(6)
!       call amrex_abort("x-y plane")
    end if

  end subroutine compute_diff_wallflux

end module cns_eb_diff_wall_module
