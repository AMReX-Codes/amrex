module cns_eb_hyp_wall_module
  use amrex_fort_module, only : rt=>amrex_real
  use amrex_error_module, only : amrex_abort
  implicit none
  private
  public :: compute_hyp_wallflux

contains

  subroutine compute_hyp_wallflux (divw, i,j,k, rho, u, v, w, p, &
       axm, axp, aym, ayp, azm, azp)
    use cns_physics_module, only : gamma
    use cns_module, only : smallp, smallr, umx, umy, umz
    use riemann_module, only : analriem
    integer, intent(in) :: i,j,k
    real(rt), intent(in) :: rho, u, v, w, p, axm, axp, aym, ayp, azm, azp
    real(rt), intent(out) :: divw(5)
    
    real(rt) :: apnorm, apnorminv, anrmx, anrmy, anrmz, un
    real(rt) :: flux(1,1,1,5)

    apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2 + (azm-azp)**2)

    if (apnorm .eq. 0.d0) then
       print *, "compute_hyp_wallflux: ", i,j,k, axm, axp, aym, ayp, azm, azp
       flush(6)
       call amrex_abort("compute_hyp_wallflux: we are in trouble.")
    end if

    apnorminv = 1.d0 / apnorm
    anrmx = (axm-axp) * apnorminv  ! pointing to the wall
    anrmy = (aym-ayp) * apnorminv
    anrmz = (azm-azp) * apnorminv

    un = u*anrmx + v*anrmy + w*anrmz

    call analriem(gamma, smallp, smallr, 1, 1, 1, 1, &
         [rho], [ un], [p], [0.d0], [0.d0], &  ! fluid
         [rho], [-un], [p], [0.d0], [0.d0], &  ! body
         flux, [1,1,1], [1,1,1], 2,3,4)

    divw = 0.d0
    divw(umx) = (axm-axp) * flux(1,1,1,2)
    divw(umy) = (aym-ayp) * flux(1,1,1,2)
    divw(umz) = (azm-azp) * flux(1,1,1,2)

  end subroutine compute_hyp_wallflux

end module cns_eb_hyp_wall_module
