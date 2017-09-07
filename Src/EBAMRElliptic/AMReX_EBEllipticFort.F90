#include "AMReX_CONSTANTS.H"

module ebefnd_module

  !     since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  !     for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains
  
  integer function identitymat(i, j)
    implicit none
    integer i, j, retval
    retval = 0
    if (i.eq.j) then
       retval = 1
    endif

    identitymat = retval
  end function identitymat

    
  subroutine ebefnd_decrinvrelcoefebco( &
       relco, relco_lo, relco_hi,  relco_nco, &
       bcoef, bcoef_lo, bcoef_hi,  bcoef_nco, &
       gridlo, gridhi, &
       beta, dx, idir) &
       bind(C, name="ebefnd_decrinvrelcoefebco")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: i, j, k, relco_nco, bcoef_nco, idir
    integer      :: relco_lo(0:2),relco_hi(0:2)
    integer      :: bcoef_lo(0:2),bcoef_hi(0:2)
    integer      :: gridlo(0:2), gridhi(0:2), ii, jj, kk
    real(c_real) :: beta, dx
    real(c_real) :: bcoef(bcoef_lo(0):bcoef_hi(0),bcoef_lo(1):bcoef_hi(1),bcoef_lo(2):bcoef_hi(2), 0:bcoef_nco-1)
    real(c_real) :: relco(relco_lo(0):relco_hi(0),relco_lo(1):relco_hi(1),relco_lo(2):relco_hi(2), 0:relco_nco-1)

    ii = identitymat(idir, 0)
    jj = identitymat(idir, 1)
    kk = identitymat(idir, 2)

    do k = gridlo(2), gridhi(2)
       do j = gridlo(1), gridhi(1)
          do i = gridlo(0), gridhi(0)

             relco(i,j,k,0) = relco(i,j,k,0) &
                  - beta*( &
                  bcoef(i+ii,j+jj,k+kk,0) + &
                  bcoef(i   ,j   ,k   ,0))/(dx*dx)

          enddo
       enddo
    enddo

    return 
  end subroutine ebefnd_decrinvrelcoefebco


  subroutine ebefnd_invertlambdaebco( &
       relco, relco_lo, relco_hi,  relco_nco, &
       gridlo, gridhi, safety) &
       bind(C, name="ebefnd_invertlambdaebco")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: i, j, k, relco_nco
    integer      :: relco_lo(0:2),relco_hi(0:2)
    integer      :: gridlo(0:2), gridhi(0:2)
    real(c_real) :: safety
    real(c_real) :: relco(relco_lo(0):relco_hi(0),relco_lo(1):relco_hi(1),relco_lo(2):relco_hi(2), 0:relco_nco-1)


    do k = gridlo(2), gridhi(2)
       do j = gridlo(1), gridhi(1)
          do i = gridlo(0), gridhi(0)

             relco(i,j,k,0) = safety/relco(i,j,k,0)

          enddo
       enddo
    enddo

    return 
  end subroutine ebefnd_invertlambdaebco

  subroutine ebefnd_applyop_ebcond_nobcs( &
       lph,   lph_lo,  lph_hi,   lph_nco, &
       phi,   phi_lo,  phi_hi,   phi_nco, &
       acoe, acoe_lo, acoe_hi,  acoe_nco, &
       bco0, bco0_lo, bco0_hi,  bco0_nco, &
       bco1, bco1_lo, bco1_hi,  bco1_nco, &
       bco2, bco2_lo, bco2_hi,  bco2_nco, &
       gridlo, gridhi, &
       dx, alpha, beta) &
       bind(C, name="ebefnd_applyop_ebcond_nobcs")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: i, j, k, lph_nco, phi_nco
    integer      :: phi_lo(0:2),phi_hi(0:2)
    integer      :: lph_lo(0:2),lph_hi(0:2)
    integer      :: acoe_nco  
    integer      :: bco0_nco  
    integer      :: bco1_nco  
    integer      :: bco2_nco  

    integer      :: acoe_lo(0:2),acoe_hi(0:2)
    integer      :: bco0_lo(0:2),bco0_hi(0:2)
    integer      :: bco1_lo(0:2),bco1_hi(0:2)
    integer      :: bco2_lo(0:2),bco2_hi(0:2)

    integer      :: gridlo(0:2), gridhi(0:2)
    real(c_real) :: xterm, yterm, zterm
    real(c_real) :: alpha, beta, dx, laplphi
    real(c_real) :: phi(phi_lo(0):phi_hi(0),phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2), 0:phi_nco-1)
    real(c_real) :: lph(lph_lo(0):lph_hi(0),lph_lo(1):lph_hi(1),lph_lo(2):lph_hi(2), 0:lph_nco-1)

    real(c_real) :: acoe(acoe_lo(0):acoe_hi(0),acoe_lo(1):acoe_hi(1),acoe_lo(2):acoe_hi(2), 0:acoe_nco-1)
    real(c_real) :: bco0(bco0_lo(0):bco0_hi(0),bco0_lo(1):bco0_hi(1),bco0_lo(2):bco0_hi(2), 0:bco0_nco-1)
    real(c_real) :: bco1(bco1_lo(0):bco1_hi(0),bco1_lo(1):bco1_hi(1),bco1_lo(2):bco1_hi(2), 0:bco1_nco-1)
    real(c_real) :: bco2(bco2_lo(0):bco2_hi(0),bco2_lo(1):bco2_hi(1),bco2_lo(2):bco2_hi(2), 0:bco2_nco-1)


    do k = gridlo(2), gridhi(2)
       do j = gridlo(1), gridhi(1)
          do i = gridlo(0), gridhi(0)

             zterm = zero
             xterm = (bco0(i+1,j  ,k  ,0)*(phi(i+1,j  ,k  ,0) - phi(i  ,j  ,k  ,0))  &
                     -bco0(i  ,j  ,k  ,0)*(phi(i  ,j  ,k  ,0) - phi(i-1,j  ,k  ,0)))/dx/dx
             yterm = (bco1(i  ,j+1,k  ,0)*(phi(i  ,j+1,k  ,0) - phi(i  ,j  ,k  ,0)) &
                     -bco1(i  ,j  ,k  ,0)*(phi(i  ,j  ,k  ,0) - phi(i  ,j-1,k  ,0)))/dx/dx
#if BL_SPACEDIM==3
             zterm = (bco2(i  ,j  ,k+1,0)*(phi(i  ,j  ,k+1,0) - phi(i  ,j  ,k  ,0)) &
                     -bco2(i  ,j  ,k  ,0)*(phi(i  ,j  ,k  ,0) - phi(i  ,j  ,k-1,0)))/dx/dx
#endif
             laplphi = xterm + yterm + zterm

             lph(i,j,k,0) = alpha*acoe(i,j,k,0)*phi(i,j,k,0) + beta*laplphi

          enddo
       enddo
    enddo

    return 
  end subroutine ebefnd_applyop_ebcond_nobcs


  subroutine ebefnd_gscolor_ebcond( &
       phi,   phi_lo,  phi_hi,   phi_nco, &
       rhs,   rhs_lo,  rhs_hi,   rhs_nco, &
       lam,   lam_lo,  lam_hi,   lam_nco, &
       acoe, acoe_lo, acoe_hi,  acoe_nco, &
       bco0, bco0_lo, bco0_hi,  bco0_nco, &
       bco1, bco1_lo, bco1_hi,  bco1_nco, &
       bco2, bco2_lo, bco2_hi,  bco2_nco, &
       gridlo, gridhi, &
       dx, alpha, beta) &
       bind(C, name="ebefnd_gscolor_ebcond")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: i, j, k, rhs_nco, phi_nco, lam_nco
    integer      :: phi_lo(0:2),phi_hi(0:2)
    integer      :: rhs_lo(0:2),rhs_hi(0:2)
    integer      :: lam_lo(0:2),lam_hi(0:2)
    integer      :: acoe_nco  
    integer      :: bco0_nco  
    integer      :: bco1_nco  
    integer      :: bco2_nco  

    integer      :: acoe_lo(0:2),acoe_hi(0:2)
    integer      :: bco0_lo(0:2),bco0_hi(0:2)
    integer      :: bco1_lo(0:2),bco1_hi(0:2)
    integer      :: bco2_lo(0:2),bco2_hi(0:2)

    integer      :: gridlo(0:2), gridhi(0:2)
    real(c_real) :: alpha, beta, dx, laplphi
    real(c_real) :: xterm, yterm, zterm, acopt, phipt, oper
    real(c_real) :: phi(phi_lo(0):phi_hi(0),phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2), 0:phi_nco-1)
    real(c_real) :: rhs(rhs_lo(0):rhs_hi(0),rhs_lo(1):rhs_hi(1),rhs_lo(2):rhs_hi(2), 0:rhs_nco-1)
    real(c_real) :: lam(lam_lo(0):lam_hi(0),lam_lo(1):lam_hi(1),lam_lo(2):lam_hi(2), 0:lam_nco-1)

    real(c_real) :: acoe(acoe_lo(0):acoe_hi(0),acoe_lo(1):acoe_hi(1),acoe_lo(2):acoe_hi(2), 0:acoe_nco-1)
    real(c_real) :: bco0(bco0_lo(0):bco0_hi(0),bco0_lo(1):bco0_hi(1),bco0_lo(2):bco0_hi(2), 0:bco0_nco-1)
    real(c_real) :: bco1(bco1_lo(0):bco1_hi(0),bco1_lo(1):bco1_hi(1),bco1_lo(2):bco1_hi(2), 0:bco1_nco-1)
    real(c_real) :: bco2(bco2_lo(0):bco2_hi(0),bco2_lo(1):bco2_hi(1),bco2_lo(2):bco2_hi(2), 0:bco2_nco-1)


    do k = gridlo(2), gridhi(2), 2
       do j = gridlo(1), gridhi(1), 2
          do i = gridlo(0), gridhi(0), 2

             zterm = zero
             xterm = (bco0(i+1,j  ,k  ,0)*(phi(i+1,j  ,k  ,0) - phi(i  ,j  ,k  ,0))  &
                     -bco0(i  ,j  ,k  ,0)*(phi(i  ,j  ,k  ,0) - phi(i-1,j  ,k  ,0)))/dx/dx
             yterm = (bco1(i  ,j+1,k  ,0)*(phi(i  ,j+1,k  ,0) - phi(i  ,j  ,k  ,0)) &
                     -bco1(i  ,j  ,k  ,0)*(phi(i  ,j  ,k  ,0) - phi(i  ,j-1,k  ,0)))/dx/dx
#if BL_SPACEDIM==3
             zterm = (bco2(i  ,j  ,k+1,0)*(phi(i  ,j  ,k+1,0) - phi(i  ,j  ,k  ,0)) &
                     -bco2(i  ,j  ,k  ,0)*(phi(i  ,j  ,k  ,0) - phi(i  ,j  ,k-1,0)))/dx/dx
#endif
             laplphi = xterm + yterm + zterm
             acopt  = acoe(i,j,k,0)
             phipt  = phi( i,j,k,0)
             oper = alpha*acopt*phipt + beta*laplphi
!             oper  = alpha*acoe(i,j,k,0)*phi(i,j,k,0) + beta*laplphi

             phi(i,j,k,0) = phi(i,j,k,0) + lam(i,j,k,0)*(rhs(i,j,k,0) - oper)

          enddo
       enddo
    enddo

    return 
  end subroutine ebefnd_gscolor_ebcond


  subroutine ebefnd_getflux_ebco( &
       flux,  flux_lo,  flux_hi,   flux_nco, &
       phi,   phi_lo,   phi_hi,    phi_nco, &
       bcoef, bcoef_lo, bcoef_hi,  bcoef_nco, &
       gridlo, gridhi, &
       dx, idir) &
       bind(C, name="ebefnd_getflux_ebco")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: phi_lo(0:2),phi_hi(0:2), phi_nco
    integer      :: i, j, k, flux_nco, bcoef_nco, idir
    integer      :: flux_lo(0:2),flux_hi(0:2)
    integer      :: bcoef_lo(0:2),bcoef_hi(0:2)
    integer      :: gridlo(0:2), gridhi(0:2), ii, jj, kk
    real(c_real) :: beta, dx
    real(c_real) :: phi(phi_lo(0):phi_hi(0),phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2), 0:phi_nco-1)
    real(c_real) :: bcoef(bcoef_lo(0):bcoef_hi(0),bcoef_lo(1):bcoef_hi(1),bcoef_lo(2):bcoef_hi(2), 0:bcoef_nco-1)
    real(c_real) :: flux(flux_lo(0):flux_hi(0),flux_lo(1):flux_hi(1),flux_lo(2):flux_hi(2), 0:flux_nco-1)

    ii = identitymat(idir, 0)
    jj = identitymat(idir, 1)
    kk = identitymat(idir, 2)

    do k = gridlo(2), gridhi(2)
       do j = gridlo(1), gridhi(1)
          do i = gridlo(0), gridhi(0)

             flux(i,j,k,0) = (beta*bcoef(i,j,k,0)/dx)*( &
                  phi(i   ,j   ,k   ,0) + &
                  phi(i-ii,j-jj,k-kk,0))

          enddo
       enddo
    enddo

    return 
  end subroutine ebefnd_getflux_ebco

end module ebefnd_module

