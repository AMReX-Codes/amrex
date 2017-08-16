module cns_eb_hyp_wall_module
  use amrex_fort_module, only : rt=>amrex_real
  use amrex_error_module, only : amrex_abort
  implicit none
  private
  public :: compute_hyp_wallflux

contains

  subroutine compute_hyp_wallflux (divw, rho, u, v, w, p, &
       axm, axp, aym, ayp, azm, azp)
    use cns_physics_module, only : gamma
    use cns_module, only : smallp, smallr, umx, umy, umz
    use riemann_module, only : analriem
    real(rt), intent(in) :: rho, u, v, w, p, axm, axp, aym, ayp, azm, azp
    real(rt), intent(out) :: divw(5)
    
    real(rt) :: apnorm, apnorminv, anrmx, anrmy, anrmz
    real(rt) :: sfluid(5), sbody(5), flux(5)

    apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2 + (azm-azp)**2)

    if (apnorm .eq. 0.d0) then
       call amrex_abort("compute_hyp_wallflux: we are in trouble.")
    end if

    apnorminv = 1.d0 / apnorm
    anrmx = (axm-axp) * apnorminv  ! pointing to the wall
    anrmy = (aym-ayp) * apnorminv
    anrmz = (azm-azp) * apnorminv

    sfluid(1) = rho
    sfluid(2) = u*anrmx + v*anrmy + w*anrmz
    sfluid(3) = p
    sfluid(4) = 0.d0
    sfluid(5) = 0.d0
    
    sbody(1) =  rho
    sbody(2) = -sfluid(2)
    sbody(3) =  p
    sbody(4) =  0.d0
    sbody(5) =  0.d0

    call analriem(gamma, sfluid, sbody, smallp, smallr, flux)

    divw = 0.d0
    divw(umx) = (axm-axp) * flux(2)
    divw(umy) = (aym-ayp) * flux(2)
    divw(umz) = (azm-azp) * flux(2)

  end subroutine compute_hyp_wallflux

end module cns_eb_hyp_wall_module
