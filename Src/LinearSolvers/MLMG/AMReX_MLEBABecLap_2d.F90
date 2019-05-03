
module amrex_mlebabeclap_2d_module

  use amrex_error_module
  use amrex_constants_module, only : zero, one, two, half, third, fourth
  use amrex_fort_module, only : amrex_real
  use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, is_single_valued_cell, &
       get_neighbor_cells_int_single
  implicit none

  private
  public :: amrex_mlebabeclap_adotx, amrex_mlebabeclap_gsrb, amrex_mlebabeclap_normalize, &
            amrex_eb_mg_interp, amrex_mlebabeclap_flux, amrex_mlebabeclap_grad, amrex_blend_beta, &
            compute_dphidn_2d, compute_dphidn_2d_ho, amrex_get_dx_eb

contains

  elemental function amrex_get_dx_eb (kappa)
    real(amrex_real), intent(in) :: kappa
    real(amrex_real) :: amrex_get_dx_eb
    amrex_get_dx_eb = max(0.3d0, (kappa*kappa-0.25d0)/(2.d0*kappa))
  end function amrex_get_dx_eb

  pure function amrex_blend_beta (kappa) result(beta)
    real(amrex_real), intent(in) :: kappa
    real(amrex_real) :: beta
#if 1
    real(amrex_real),parameter :: blend_kappa = -1.d-4
    if (kappa .lt. blend_kappa) then
       beta = zero
    else
       beta = one
    end if
#else
    real(amrex_real),parameter :: blend_kappa = half
    real(amrex_real),parameter :: kapinv = one/(one-blend_kappa)
    beta = kapinv*(kappa-blend_kappa)
    beta = min(one,max(zero,beta))
#endif
  end function amrex_blend_beta

  subroutine amrex_mlebabeclap_adotx(lo, hi, y, ylo, yhi, x, xlo, xhi, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, ccm, cmlo, cmhi, flag, flo, fhi, vfrc, vlo, vhi, &
       apx, axlo, axhi, apy, aylo, ayhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
       ba, balo, bahi, bc, bclo, bchi, beb, elo, ehi, &
       is_dirichlet, is_ho_dirichlet, phieb, plo, phi, is_inhomog, dxinv, alpha, beta) &
       bind(c,name='amrex_mlebabeclap_adotx')

    integer, dimension(2), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, alo, ahi, bxlo, bxhi, bylo, byhi, &
         cmlo, cmhi, flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, cxlo, cxhi, cylo, cyhi, balo, bahi, &
         bclo, bchi, elo, ehi, plo, phi
    real(amrex_real), intent(in) :: dxinv(2)
    integer         , value, intent(in) :: is_dirichlet, is_ho_dirichlet, is_inhomog
    real(amrex_real), value, intent(in) :: alpha, beta
    real(amrex_real), intent(inout) ::    y( ylo(1): yhi(1), ylo(2): yhi(2))
    real(amrex_real), intent(in   ) ::    x( xlo(1): xhi(1), xlo(2): xhi(2))
    real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2))
    real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2))
    integer         , intent(in   ) ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2))
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2))
    real(amrex_real), intent(in   ) :: vfrc( vlo(1): vhi(1), vlo(2): vhi(2))
    real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2))
    real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2))
    real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2))
    real(amrex_real), intent(in   ) ::   ba(balo(1):bahi(1),balo(2):bahi(2))
    real(amrex_real), intent(in   ) ::   bc(bclo(1):bchi(1),bclo(2):bchi(2),2)
    real(amrex_real), intent(in   ) ::  beb( elo(1): ehi(1), elo(2): ehi(2))
    real(amrex_real), intent(in   ) ::phieb( plo(1): phi(1), plo(2): phi(2))
    integer :: i,j, ii, jj
    real(amrex_real) :: dhx, dhy, fxm, fxp, fym, fyp, fracx, fracy
    real(amrex_real) :: feb, phib, anrmx, anrmy, anorm, anorminv
    logical :: is_inhomogeneous

    real(amrex_real) :: dphidn

    is_inhomogeneous = is_inhomog .ne. 0

    dhx = beta*dxinv(1)*dxinv(1)
    dhy = beta*dxinv(2)*dxinv(2)
    
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (is_covered_cell(flag(i,j))) then
             y(i,j) = zero
          else if (is_regular_cell(flag(i,j))) then
             y(i,j) = alpha*a(i,j)*x(i,j) &
                  - dhx * (bX(i+1,j)*(x(i+1,j) - x(i  ,j))  &
                  &      - bX(i  ,j)*(x(i  ,j) - x(i-1,j))) &
                  - dhy * (bY(i,j+1)*(x(i,j+1) - x(i,j  ))  &
                  &      - bY(i,j  )*(x(i,j  ) - x(i,j-1)))
          else
             fxm = bX(i,j)*(x(i,j)-x(i-1,j))
             if (apx(i,j).ne.zero .and. apx(i,j).ne.one) then
                jj = j + int(sign(one,fcx(i,j)))
                fracy = abs(fcx(i,j))*real(ior(ccm(i-1,jj),ccm(i,jj)),amrex_real)
                fxm = (one-fracy)*fxm + fracy*bX(i,jj)*(x(i,jj)-x(i-1,jj))
             end if

             fxp = bX(i+1,j)*(x(i+1,j)-x(i,j))
             if (apx(i+1,j).ne.zero .and. apx(i+1,j).ne.one) then
                jj = j + int(sign(one,fcx(i+1,j)))
                fracy = abs(fcx(i+1,j))*real(ior(ccm(i,jj),ccm(i+1,jj)),amrex_real)
                fxp = (one-fracy)*fxp + fracy*bX(i+1,jj)*(x(i+1,jj)-x(i,jj))
             end if

             fym = bY(i,j)*(x(i,j)-x(i,j-1))
             if (apy(i,j).ne.zero .and. apy(i,j).ne.one) then
                ii = i + int(sign(one,fcy(i,j)))
                fracx = abs(fcy(i,j))*real(ior(ccm(ii,j-1),ccm(ii,j)),amrex_real)
                fym = (one-fracx)*fym + fracx*bY(ii,j)*(x(ii,j)-x(ii,j-1))
             end if

             fyp = bY(i,j+1)*(x(i,j+1)-x(i,j))
             if (apy(i,j+1).ne.zero .and. apy(i,j+1).ne.one) then
                ii = i + int(sign(one,fcy(i,j+1)))
                fracx = abs(fcy(i,j+1))*real(ior(ccm(ii,j),ccm(ii,j+1)),amrex_real)
                fyp = (one-fracx)*fyp + fracx*bY(ii,j+1)*(x(ii,j+1)-x(ii,j))
             end if

             if (is_ho_dirichlet .ne. 0 .or. is_dirichlet .ne. 0) then
                anorm = sqrt((apx(i,j)-apx(i+1,j))**2 + (apy(i,j)-apy(i,j+1))**2)
                anorminv = one/anorm
                anrmx = (apx(i,j)-apx(i+1,j)) * anorminv
                anrmy = (apy(i,j)-apy(i,j+1)) * anorminv

                if (is_inhomogeneous) then
                   phib = phieb(i,j)
                else
                   phib = zero
                end if

                if (is_ho_dirichlet .ne. 0) then
                   call compute_dphidn_2d_ho(dphidn, dxinv, i, j,  &
                                             x,    xlo,  xhi, &
                                             flag, flo,  fhi, &
                                             bc(i,j,:),  phib, &
                                             anrmx, anrmy)
                else if (is_dirichlet .ne. 0) then
                   call compute_dphidn_2d(dphidn, dxinv, i, j,  &
                                          x,    xlo,  xhi, &
                                          flag, flo,  fhi, &
                                          bc(i,j,:),  phib, &
                                          anrmx, anrmy, vfrc(i,j))
                end if
                feb = dphidn * ba(i,j) * beb(i,j)
             else
                feb = zero
             end if

             y(i,j) = alpha*a(i,j)*x(i,j) + (one/vfrc(i,j)) * &
                  (dhx*(apx(i,j)*fxm-apx(i+1,j)*fxp) + dhy*(apy(i,j)*fym-apy(i,j+1)*fyp) &
                  -dhx*feb)
          end if
       end do
    end do
  end subroutine amrex_mlebabeclap_adotx

  subroutine amrex_mlebabeclap_gsrb(lo, hi, phi, hlo, hhi, rhs, rlo, rhi, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, &
       ccm, cmlo, cmhi, &
       m0, m0lo, m0hi, m2, m2lo, m2hi, &
       m1, m1lo, m1hi, m3, m3lo, m3hi, &
       f0, f0lo, f0hi, f2, f2lo, f2hi, &
       f1, f1lo, f1hi, f3, f3lo, f3hi, &
       flag, flo, fhi, vfrc, vlo, vhi, &
       apx, axlo, axhi, apy, aylo, ayhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
       ba, balo, bahi, bc, bclo, bchi, beb, elo, ehi, &
       is_dirichlet, is_ho_dirichlet, &
       dxinv, alpha, beta, redblack) &
       bind(c,name='amrex_mlebabeclap_gsrb')

    integer, dimension(2), intent(in) :: lo, hi, hlo, hhi, rlo, rhi, alo, ahi, bxlo, bxhi, bylo, byhi, &
         cmlo, cmhi, m0lo, m0hi, m1lo, m1hi, m2lo, m2hi, m3lo, m3hi, &
         f0lo, f0hi, f1lo, f1hi, f2lo, f2hi, f3lo, f3hi, &
         flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, cxlo, cxhi, cylo, cyhi, &
         balo, bahi, bclo, bchi, elo, ehi

    real(amrex_real), intent(in) :: dxinv(2)
    integer         , value, intent(in) :: is_dirichlet, is_ho_dirichlet
    real(amrex_real), value, intent(in) :: alpha, beta
    integer, value, intent(in) :: redblack
    real(amrex_real), intent(inout) ::  phi( hlo(1): hhi(1), hlo(2): hhi(2))
    real(amrex_real), intent(in   ) ::  rhs( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2))
    real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2))
    integer         , intent(in   ) ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2))
    integer         , intent(in   ) ::   m0(m0lo(1):m0hi(1),m0lo(2):m0hi(2))
    integer         , intent(in   ) ::   m1(m1lo(1):m1hi(1),m1lo(2):m1hi(2))
    integer         , intent(in   ) ::   m2(m2lo(1):m2hi(1),m2lo(2):m2hi(2))
    integer         , intent(in   ) ::   m3(m3lo(1):m3hi(1),m3lo(2):m3hi(2))
    real(amrex_real), intent(in   ) ::   f0(f0lo(1):f0hi(1),f0lo(2):f0hi(2))
    real(amrex_real), intent(in   ) ::   f1(f1lo(1):f1hi(1),f1lo(2):f1hi(2))
    real(amrex_real), intent(in   ) ::   f2(f2lo(1):f2hi(1),f2lo(2):f2hi(2))
    real(amrex_real), intent(in   ) ::   f3(f3lo(1):f3hi(1),f3lo(2):f3hi(2))
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2))
    real(amrex_real), intent(in   ) :: vfrc( vlo(1): vhi(1), vlo(2): vhi(2))
    real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2))
    real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2))
    real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2))
    real(amrex_real), intent(in   ) ::   ba(balo(1):bahi(1),balo(2):bahi(2))
    real(amrex_real), intent(in   ) ::   bc(bclo(1):bchi(1),bclo(2):bchi(2),2)
    real(amrex_real), intent(in   ) ::  beb( elo(1): ehi(1), elo(2): ehi(2))

    integer :: i,j,ioff,ii,jj
    real(amrex_real) :: cf0, cf1, cf2, cf3, delta, gamma, rho, res, vfrcinv
    real(amrex_real) :: dhx, dhy, fxm, fxp, fym, fyp, fracx, fracy
    real(amrex_real) :: sxm, sxp, sym, syp, oxm, oxp, oym, oyp
    real(amrex_real) :: feb, phib, phig, phig1, phig2, gx, gy, anrmx, anrmy, anorm, anorminv, sx, sy
    real(amrex_real) :: feb_gamma, phig_gamma, phig1_gamma
    real(amrex_real) :: bctx, bcty, bsxinv, bsyinv
    real(amrex_real) :: w1, w2, dg, dx_eb
    real(amrex_real), dimension(-1:0,-1:0) :: c_0, c_x, c_y, c_xy
    real(amrex_real), parameter :: omega = 1._amrex_real

    real(amrex_real) :: dphidn

    dhx = beta*dxinv(1)*dxinv(1)
    dhy = beta*dxinv(2)*dxinv(2)

    do j = lo(2), hi(2)
       ioff = mod(lo(1)+j+redblack,2)
       do i = lo(1)+ioff, hi(1), 2

          if (is_covered_cell(flag(i,j))) then
             phi(i,j) = zero
          else
             cf0 = merge(f0(lo(1),j), 0.0D0, &
                  (i .eq. lo(1)) .and. (m0(lo(1)-1,j).gt.0))
             cf1 = merge(f1(i,lo(2)), 0.0D0, &
                  (j .eq. lo(2)) .and. (m1(i,lo(2)-1).gt.0))
             cf2 = merge(f2(hi(1),j), 0.0D0, &
                  (i .eq. hi(1)) .and. (m2(hi(1)+1,j).gt.0))
             cf3 = merge(f3(i,hi(2)), 0.0D0, &
                  (j .eq. hi(2)) .and. (m3(i,hi(2)+1).gt.0))
             
             if (is_regular_cell(flag(i,j))) then
                
                gamma = alpha*a(i,j) &
                     + dhx * (bX(i+1,j) + bX(i,j)) &
                     + dhy * (bY(i,j+1) + bY(i,j))
                
                rho =  dhx * (bX(i+1,j)*phi(i+1,j) + bX(i,j)*phi(i-1,j)) &
                     + dhy * (bY(i,j+1)*phi(i,j+1) + bY(i,j)*phi(i,j-1))

                delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2) &
                     +  dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)
             
             else

                fxm = -bX(i,j)*phi(i-1,j)
                oxm = -bX(i,j)*cf0
                sxm =  bX(i,j)
                if (apx(i,j).ne.zero .and. apx(i,j).ne.one) then
                   jj = j + int(sign(one,fcx(i,j)))
                   fracy = abs(fcx(i,j))*real(ior(ccm(i-1,jj),ccm(i,jj)),amrex_real)
                   fxm = (one-fracy)*fxm + fracy*bX(i,jj)*(phi(i,jj)-phi(i-1,jj))
                   ! oxm = (one-fracy)*oxm
                   oxm = zero
                   sxm = (one-fracy)*sxm
                end if
                
                fxp =  bX(i+1,j)*phi(i+1,j)
                oxp =  bX(i+1,j)*cf2
                sxp = -bX(i+1,j)
                if (apx(i+1,j).ne.zero .and. apx(i+1,j).ne.one) then
                   jj = j + int(sign(one,fcx(i+1,j)))
                   fracy = abs(fcx(i+1,j))*real(ior(ccm(i,jj),ccm(i+1,jj)),amrex_real)
                   fxp = (one-fracy)*fxp + fracy*bX(i+1,jj)*(phi(i+1,jj)-phi(i,jj))
                   ! oxp = (one-fracy)*oxp
                   oxp = zero
                   sxp = (one-fracy)*sxp
                end if
                
                fym = -bY(i,j)*phi(i,j-1)
                oym = -bY(i,j)*cf1
                sym =  bY(i,j)
                if (apy(i,j).ne.zero .and. apy(i,j).ne.one) then
                   ii = i + int(sign(one,fcy(i,j)))
                   fracx = abs(fcy(i,j))*real(ior(ccm(ii,j-1),ccm(ii,j)),amrex_real)
                   fym = (one-fracx)*fym + fracx*bY(ii,j)*(phi(ii,j)-phi(ii,j-1))
                   ! oym = (one-fracx)*oym
                   oym = zero
                   sym = (one-fracx)*sym
                end if
                
                fyp =  bY(i,j+1)*phi(i,j+1)
                oyp =  bY(i,j+1)*cf3
                syp = -bY(i,j+1)
                if (apy(i,j+1).ne.zero .and. apy(i,j+1).ne.one) then
                   ii = i + int(sign(one,fcy(i,j+1)))
                   fracx = abs(fcy(i,j+1))*real(ior(ccm(ii,j),ccm(ii,j+1)),amrex_real)
                   fyp = (one-fracx)*fyp + fracx*bY(ii,j+1)*(phi(ii,j+1)-phi(ii,j))
                   ! oyp = (one-fracx)*fyp
                   oyp = zero
                   syp = (one-fracx)*syp
                end if

                vfrcinv = (one/vfrc(i,j))
                gamma = alpha*a(i,j) + vfrcinv * &
                     (dhx*(apx(i,j)*sxm-apx(i+1,j)*sxp) + dhy*(apy(i,j)*sym-apy(i,j+1)*syp))
                rho = -vfrcinv * &
                     (dhx*(apx(i,j)*fxm-apx(i+1,j)*fxp) + dhy*(apy(i,j)*fym-apy(i,j+1)*fyp))

                delta = -vfrcinv * &
                     (dhx*(apx(i,j)*oxm-apx(i+1,j)*oxp) + dhy*(apy(i,j)*oym-apy(i,j+1)*oyp))

                if (is_ho_dirichlet .ne. 0 .or. is_dirichlet .ne. 0) then

                   anorm = sqrt((apx(i,j)-apx(i+1,j))**2 + (apy(i,j)-apy(i,j+1))**2)
                   anorminv = one/anorm
                   anrmx = (apx(i,j)-apx(i+1,j)) * anorminv
                   anrmy = (apy(i,j)-apy(i,j+1)) * anorminv

                   ! In gsrb we are always in residual-correction form so phib = 0
                   phib = zero

                end if

                if (is_ho_dirichlet .ne. 0) then
                    call compute_dphidn_2d_ho(dphidn, dxinv, i, j,  &
                                              phi,  hlo,  hhi, &
                                              flag, flo,  fhi, &
                                              bc(i,j,:),  phib, &
                                              anrmx, anrmy)

                   ! We should modify these but haven't done it yet
                   ! feb_gamma = -phig_gamma * (ba(i,j)*beb(i,j)/dg)
                   ! gamma = gamma + vfrcinv*(-dhx)*feb_gamma

                   feb = dphidn * ba(i,j) * beb(i,j)
                   rho = rho - vfrcinv*(-dhx)*feb

                else if (is_dirichlet .ne. 0) then

                   bctx = bc(i,j,1)
                   bcty = bc(i,j,2)

                   dx_eb = amrex_get_dx_eb(vfrc(i,j))

                   if (abs(anrmx) .gt. abs(anrmy)) then
                      dg = dx_eb / abs(anrmx)
                      gx = bctx - dg*anrmx
                      gy = bcty - dg*anrmy
                      sx =  sign(one,anrmx)
                      sy =  sign(one,anrmy)
                      ! sy = -sign(one,gy)
                   else
                      dg = dx_eb / abs(anrmy)
                      gx = bctx - dg*anrmx
                      gy = bcty - dg*anrmy
                      ! sx = -sign(one,gx)
                      sx =  sign(one,anrmx)
                      sy =  sign(one,anrmy)
                   end if
                   ii = i - int(sx)
                   jj = j - int(sy)
                   
                   w1 = amrex_blend_beta(vfrc(i,j))
                   w2 = one-w1

                   if (w1.eq.zero) then
                      phig1_gamma = zero
                      phig1 = zero
                   else
                      phig1_gamma = (one + gx*sx + gy*sy + gx*gy*sx*sy)
                      phig1 = (    - gx*sx         - gx*gy*sx*sy) * phi(ii,j) &
                           +  (            - gy*sy - gx*gy*sx*sy) * phi(i,jj) &
                           +  (                    + gx*gy*sx*sy) * phi(ii,jj)
                   end if

                   if (w2.eq.zero) then
                      phig2 = zero
                   else
                      bsxinv = one/(bctx+sx)
                      bsyinv = one/(bcty+sy)
                   
                      ! c_0(0,0) = sx*sy*bsxinv*bsyinv
                      c_0(-1,0) = bctx*bsxinv
                      c_0(0,-1) = bcty*bsyinv
                      c_0(-1,-1) = -bctx*bcty*bsxinv*bsyinv
                   
                      ! c_x(0,0) = sy*bsxinv*bsyinv
                      c_x(-1,0) = -bsxinv
                      c_x(0,-1) = sx*bcty*bsyinv
                      c_x(-1,-1) = -sx*bctx*bcty*bsxinv*bsyinv
                      
                      ! c_y(0,0) = sx*bsxinv*bsyinv
                      c_y(-1,0) = sy*bctx*bsxinv
                      c_y(0,-1) = -bsyinv
                      c_y(-1,-1) = -sy*bctx*bcty*bsxinv*bsyinv
                      
                      ! c_xy(0,0) = bsxinv*bsyinv
                      c_xy(-1,0) = -sy*bsxinv
                      c_xy(0,-1) = -sx*bsyinv
                      c_xy(-1,-1) = (one+sx*bctx+sy*bcty)*bsxinv*bsyinv
                      
                      phig2 = (c_0(-1, 0) + gx*c_x(-1, 0) + gy*c_y(-1, 0) + gx*gy*c_xy(-1, 0))*phi(ii,j) &
                           +  (c_0( 0,-1) + gx*c_x( 0,-1) + gy*c_y( 0,-1) + gx*gy*c_xy( 0,-1))*phi(i,jj) &
                           +  (c_0(-1,-1) + gx*c_x(-1,-1) + gy*c_y(-1,-1) + gx*gy*c_xy(-1,-1))*phi(ii,jj)
                   end if

                   phig_gamma = w1*phig1_gamma
                   phig = w1*phig1 + w2*phig2

                   dphidn =  (    -phig)/dg

                   feb = dphidn * ba(i,j) * beb(i,j)
                   rho = rho - vfrcinv*(-dhx)*feb

                   feb_gamma = -phig_gamma/dg * ba(i,j) * beb(i,j)
                   gamma = gamma + vfrcinv*(-dhx)*feb_gamma
                end if

             end if

             res = rhs(i,j) - (gamma*phi(i,j) - rho)
             phi(i,j) = phi(i,j) + omega*res/(gamma-delta)

          end if
       end do
    end do

  end subroutine amrex_mlebabeclap_gsrb


  subroutine amrex_mlebabeclap_normalize (lo, hi, x, xlo, xhi, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, ccm, cmlo, cmhi, flag, flo, fhi, vfrc, vlo, vhi, &
       apx, axlo, axhi, apy, aylo, ayhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
       ba, balo, bahi, bc, bclo, bchi, beb, elo, ehi, is_dirichlet, is_ho_dirichlet, dxinv, alpha, beta) &
       bind(c,name='amrex_mlebabeclap_normalize')

    integer, dimension(2), intent(in) :: lo, hi, xlo, xhi, alo, ahi, bxlo, bxhi, bylo, byhi, &
         cmlo, cmhi, flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, cxlo, cxhi, cylo, cyhi, &
         balo, bahi, bclo, bchi, elo, ehi
    real(amrex_real), intent(in) :: dxinv(2)
    integer         , value, intent(in) :: is_dirichlet, is_ho_dirichlet
    real(amrex_real), value, intent(in) :: alpha, beta
    real(amrex_real), intent(inout) ::    x( xlo(1): xhi(1), xlo(2): xhi(2))
    real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2))
    real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2))
    integer         , intent(in   ) ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2))
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2))
    real(amrex_real), intent(in   ) :: vfrc( vlo(1): vhi(1), vlo(2): vhi(2))
    real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2))
    real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2))
    real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2))
    real(amrex_real), intent(in   ) ::   ba(balo(1):bahi(1),balo(2):bahi(2))
    real(amrex_real), intent(in   ) ::   bc(bclo(1):bchi(1),bclo(2):bchi(2),2)
    real(amrex_real), intent(in   ) ::  beb( elo(1): ehi(1), elo(2): ehi(2))

    integer :: i,j,ii,jj
    real(amrex_real) :: dhx, dhy, sxm, sxp, sym, syp, gamma, fracx, fracy, vfrcinv
    real(amrex_real) :: gx, gy, anrmx, anrmy, anorm, anorminv, sx, sy
    real(amrex_real) :: feb_gamma, phig_gamma, phig1_gamma
    real(amrex_real) :: bctx, bcty
    real(amrex_real) :: w1, w2, dg, dx_eb

    dhx = beta*dxinv(1)*dxinv(1)
    dhy = beta*dxinv(2)*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (is_regular_cell(flag(i,j))) then
             x(i,j) = x(i,j) / (alpha*a(i,j) + dhx*(bX(i,j)+bX(i+1,j)) &
                  &                          + dhy*(bY(i,j)+bY(i,j+1)))
          else if (is_single_valued_cell(flag(i,j))) then

             sxm =  bX(i,j)
             if (apx(i,j).ne.zero .and. apx(i,j).ne.one) then
                jj = j + int(sign(one,fcx(i,j)))
                fracy = abs(fcx(i,j))*real(ior(ccm(i-1,jj),ccm(i,jj)),amrex_real)
                sxm = (one-fracy)*sxm
             end if
                
             sxp = -bX(i+1,j)
             if (apx(i+1,j).ne.zero .and. apx(i+1,j).ne.one) then
                jj = j + int(sign(one,fcx(i+1,j)))
                fracy = abs(fcx(i+1,j))*real(ior(ccm(i,jj),ccm(i+1,jj)),amrex_real)
                sxp = (one-fracy)*sxp
             end if
                
             sym =  bY(i,j)
             if (apy(i,j).ne.zero .and. apy(i,j).ne.one) then
                ii = i + int(sign(one,fcy(i,j)))
                fracx = abs(fcy(i,j))*real(ior(ccm(ii,j-1),ccm(ii,j)),amrex_real)
                sym = (one-fracx)*sym
             end if
                
             syp = -bY(i,j+1)
             if (apy(i,j+1).ne.zero .and. apy(i,j+1).ne.one) then
                ii = i + int(sign(one,fcy(i,j+1)))
                fracx = abs(fcy(i,j+1))*real(ior(ccm(ii,j),ccm(ii,j+1)),amrex_real)
                syp = (one-fracx)*syp
             end if

             vfrcinv = one/vfrc(i,j)
             gamma = alpha*a(i,j) + vfrcinv * &
                  (dhx*(apx(i,j)*sxm-apx(i+1,j)*sxp) + dhy*(apy(i,j)*sym-apy(i,j+1)*syp))

             if (is_dirichlet .ne. 0 .or. is_ho_dirichlet .ne. 0) then
                anorm = sqrt((apx(i,j)-apx(i+1,j))**2 + (apy(i,j)-apy(i,j+1))**2)
                anorminv = one/anorm
                anrmx = (apx(i,j)-apx(i+1,j)) * anorminv
                anrmy = (apy(i,j)-apy(i,j+1)) * anorminv
                bctx = bc(i,j,1)
                bcty = bc(i,j,2)

                dx_eb = amrex_get_dx_eb(vfrc(i,j))

                if (abs(anrmx) .gt. abs(anrmy)) then
                   dg = dx_eb / abs(anrmx)
                   gx = bctx - dg*anrmx
                   gy = bcty - dg*anrmy
                   sx =  sign(one,anrmx)
                   sy =  sign(one,anrmy)
                else
                   dg = dx_eb / abs(anrmy)
                   gx = bctx - dg*anrmx
                   gy = bcty - dg*anrmy
                   sx =  sign(one,anrmx)
                   sy =  sign(one,anrmy)
                end if
                ii = i - int(sx)
                jj = j - int(sy)
                
                w1 = amrex_blend_beta(vfrc(i,j))
                w2 = one-w1
                
                if (w1.eq.zero) then
                   phig1_gamma = zero
                else
                   phig1_gamma = (one + gx*sx + gy*sy + gx*gy*sx*sy)
                end if
                
                phig_gamma = w1*phig1_gamma
                feb_gamma = -phig_gamma * (ba(i,j) * beb(i,j) / dg)
                
                gamma = gamma + vfrcinv*(-dhx)*feb_gamma
             end if

             x(i,j) = x(i,j) / gamma
          end if
       end do
    end do

  end subroutine amrex_mlebabeclap_normalize

  subroutine amrex_eb_mg_interp (lo, hi, fine, flo, fhi, crse, clo, chi, flag, glo, ghi, ncomp) &
       bind(c,name='amrex_eb_mg_interp')
    integer, dimension(2), intent(in) :: lo, hi, flo, fhi, clo, chi, glo, ghi
    integer, intent(in) :: ncomp
    real(amrex_real), intent(inout) :: fine(flo(1):fhi(1),flo(2):fhi(2),ncomp)
    real(amrex_real), intent(in   ) :: crse(clo(1):chi(1),clo(2):chi(2),ncomp)
    integer         , intent(in   ) :: flag(glo(1):ghi(1),glo(2):ghi(2))

    integer :: i,j,ii,jj,n

    do n = 1, ncomp
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ii = 2*i
             jj = 2*j
             if (.not.is_covered_cell(flag(ii,jj))) then
                fine(ii,jj,n) = fine(ii,jj,n) + crse(i,j,n)
             end if

             ii = 2*i+1
             jj = 2*j
             if (.not.is_covered_cell(flag(ii,jj))) then
                fine(ii,jj,n) = fine(ii,jj,n) + crse(i,j,n)
             end if

             ii = 2*i
             jj = 2*j+1
             if (.not.is_covered_cell(flag(ii,jj))) then
                fine(ii,jj,n) = fine(ii,jj,n) + crse(i,j,n)
             end if

             ii = 2*i+1
             jj = 2*j+1
             if (.not.is_covered_cell(flag(ii,jj))) then
                fine(ii,jj,n) = fine(ii,jj,n) + crse(i,j,n)
             end if

          end do
       end do
    end do

  end subroutine amrex_eb_mg_interp

  subroutine amrex_mlebabeclap_flux(lo, hi, fx, fxlo, fxhi, fy, fylo, fyhi,  apx, axlo, axhi, & 
                                    apy, aylo, ayhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
                                    sol, slo, shi, bx, bxlo, bxhi, by, bylo, byhi,&
                                    ccm, cmlo, cmhi, flag, glo, ghi, &
                                    dxinv, beta, face_only) &
                                    bind(c, name='amrex_mlebabeclap_flux')
    integer, dimension(2), intent(in)   :: lo, hi, fxlo, fxhi, fylo, fyhi, axlo, axhi, aylo, ayhi, glo, ghi
    integer, dimension(2), intent(in)   :: cxlo, cxhi, cylo, cyhi, slo, shi, bxlo, bxhi, bylo, byhi, cmlo, cmhi

    integer,   value, intent(in   )     :: face_only
    real(amrex_real), value, intent(in) :: beta
    real(amrex_real), intent(in   )     :: dxinv(2) 
    real(amrex_real), intent(inout)     :: fx  (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(amrex_real), intent(inout)     :: fy  (fylo(1):fyhi(1),fylo(2):fyhi(2)) 
    real(amrex_real), intent(in   )     :: apx (axlo(1):axhi(1),axlo(2):axhi(2)) 
    real(amrex_real), intent(in   )     :: apy (aylo(1):ayhi(1),aylo(2):ayhi(2)) 
    real(amrex_real), intent(in   )     :: fcx (cxlo(1):cxhi(1),cxlo(2):cxhi(2))
    real(amrex_real), intent(in   )     :: fcy (cylo(1):cyhi(1),cylo(2):cyhi(2))
    real(amrex_real), intent(in   )     :: sol ( slo(1): shi(1), slo(2): shi(2))
    real(amrex_real), intent(in   )     :: bx  (bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(amrex_real), intent(in   )     :: by  (bylo(1):byhi(1),bylo(2):byhi(2)) 
    integer         , intent(in   )     :: flag( glo(1): ghi(1), glo(2): ghi(2))
    integer         , intent(in   )     ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2))

    integer :: i,j, ii, jj, istride, jstride
    real(amrex_real) :: dhx, dhy, fxm, fym, fracx, fracy

    dhx = beta*dxinv(1)
    dhy = beta*dxinv(2)

    if (face_only .eq. 1) then
       istride = hi(1)+1-lo(1)
       jstride = hi(2)+1-lo(2)
    else
       istride = 1
       jstride = 1
    end if

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1, istride
          if (apx(i,j) .eq. zero) then
             fx(i,j) = zero
          else if (is_regular_cell(flag(i,j)) .or. apx(i,j).eq.one) then
             fx(i,j) = -dhx*bx(i,j)*(sol(i,j) - sol(i-1,j))
          else
             fxm = bX(i,j)*(sol(i,j)-sol(i-1,j))
             jj = j + int(sign(one,fcx(i,j)))
             fracy = abs(fcx(i,j))*real(ior(ccm(i-1,jj),ccm(i,jj)),amrex_real)
             fxm = (one-fracy)*fxm + fracy*bX(i,jj)*(sol(i,jj)-sol(i-1,jj))
             fx(i,j) = -fxm*dhx
          end if
       end do
    end do

    do    j = lo(2), hi(2)+1, jstride
       do i = lo(1), hi(1)
          if (apy(i,j) .eq. zero) then
             fy(i,j) = zero
          else if (is_regular_cell(flag(i,j)) .or. apy(i,j).eq.one) then
             fy(i,j) = -dhy*by(i,j)*(sol(i,j) - sol(i,j-1))
          else
             fym = bY(i,j)*(sol(i,j)-sol(i,j-1))
             ii = i + int(sign(one,fcy(i,j)))
             fracx = abs(fcy(i,j))*real(ior(ccm(ii,j-1),ccm(ii,j)),amrex_real)
             fym = (one-fracx)*fym + fracx*bY(ii,j)*(sol(ii,j)-sol(ii,j-1))
             fy(i,j) = -fym*dhy
          end if
       end do
    end do
  end subroutine amrex_mlebabeclap_flux

  subroutine amrex_mlebabeclap_grad(xlo, xhi, ylo, yhi, sol, slo, shi, gx, gxlo, gxhi, & 
                                    gy, gylo, gyhi, apx, axlo, axhi, apy, aylo, ayhi,    &
                                    fcx, cxlo, cxhi, fcy, cylo, cyhi, &
                                    ccm, cmlo, cmhi, flag, glo, ghi, dxinv) &
                                    bind(c, name='amrex_mlebabeclap_grad')
    integer, dimension(2), intent(in)   :: xlo, xhi, gxlo, gxhi, gylo, gyhi, axlo, axhi, aylo, ayhi, glo, ghi 
    integer, dimension(2), intent(in)   :: ylo, yhi, cxlo, cxhi, cylo, cyhi, slo, shi, cmlo, cmhi

    real(amrex_real), intent(in   )     :: dxinv(2) 
    real(amrex_real), intent(inout)     :: gx  (gxlo(1):gxhi(1),gxlo(2):gxhi(2))
    real(amrex_real), intent(inout)     :: gy  (gylo(1):gyhi(1),gylo(2):gyhi(2)) 
    real(amrex_real), intent(in   )     :: apx (axlo(1):axhi(1),axlo(2):axhi(2)) 
    real(amrex_real), intent(in   )     :: apy (aylo(1):ayhi(1),aylo(2):ayhi(2)) 
    real(amrex_real), intent(in   )     :: fcx (cxlo(1):cxhi(1),cxlo(2):cxhi(2))
    real(amrex_real), intent(in   )     :: fcy (cylo(1):cyhi(1),cylo(2):cyhi(2))
    real(amrex_real), intent(in   )     :: sol ( slo(1): shi(1), slo(2): shi(2))
    integer         , intent(in   )     ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2))
    integer         , intent(in   )     :: flag( glo(1): ghi(1), glo(2): ghi(2))

    integer :: i,j, ii, jj
    real(amrex_real) :: dhx, dhy, fxm, fym, fracx, fracy

    dhx = dxinv(1)
    dhy = dxinv(2)
    do    j = xlo(2), xhi(2)
       do i = xlo(1), xhi(1)
          if (apx(i,j) .eq. zero) then
             gx(i,j) = zero
          else if (is_regular_cell(flag(i,j)) .or. apx(i,j).eq.one) then
             gx(i,j) = dhx*(sol(i,j) - sol(i-1,j))
          else
             fxm = (sol(i,j)-sol(i-1,j))
             jj = j + int(sign(one,fcx(i,j)))
             fracy = abs(fcx(i,j))*real(ior(ccm(i-1,jj),ccm(i,jj)),amrex_real)
             fxm = (one-fracy)*fxm + fracy*(sol(i,jj)-sol(i-1,jj))
             gx(i,j) = fxm*dhx
          end if
       end do
    end do

    do    j = ylo(2), yhi(2)
       do i = ylo(1), yhi(1)
          if (apy(i,j) .eq. zero) then
             gy(i,j) = zero
          else if (is_regular_cell(flag(i,j)) .or. apy(i,j).eq.one) then
             gy(i,j) = dhy*(sol(i,j) - sol(i,j-1))
          else
             fym = (sol(i,j)-sol(i,j-1))
             ii = i + int(sign(one,fcy(i,j)))
             fracx = abs(fcy(i,j))*real(ior(ccm(ii,j-1),ccm(ii,j)),amrex_real)
             fym = (one-fracx)*fym + fracx*(sol(ii,j)-sol(ii,j-1))
             gy(i,j) = fym*dhy
          end if
       end do
    end do
  end subroutine amrex_mlebabeclap_grad

  subroutine compute_dphidn_2d (dphidn, dxinv, i, j,  &
        phi , p_lo, p_hi,     &
        flag,  flo,  fhi,     &
        bct, phib, anrmx, anrmy, vf)

       ! Cell indices 
       integer, intent(in   ) :: i, j
 
       ! Grid spacing
       real(amrex_real),       intent(in   ) :: dxinv(2)
 
       integer, intent(in   ) :: p_lo(2), p_hi(2)
       integer, intent(in   ) ::  flo(2),  fhi(2)
       integer, intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2) )
 
       ! Arrays
       real(amrex_real),  intent(in   ) ::                            &
            & phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2))
 
       real(amrex_real),  intent(in   ) :: bct(2), phib
       real(amrex_real),  intent(in   ) :: anrmx, anrmy, vf
 
       real(amrex_real),        intent(  out) :: dphidn
 
       real(amrex_real) :: bctx, bcty, bsxinv, bsyinv
       real(amrex_real) :: w1, w2, dg, dx_eb
       real(amrex_real), dimension(-1:0,-1:0) :: c_0, c_x, c_y, c_xy
       real(amrex_real) :: phig, phig1, phig2, gx, gy, sx, sy
       integer          :: ii, jj

       bctx = bct(1)
       bcty = bct(2)

       dx_eb = amrex_get_dx_eb(vf)

       if (abs(anrmx) .gt. abs(anrmy)) then
          dg = dx_eb / abs(anrmx)
          gx = bctx - dg*anrmx
          gy = bcty - dg*anrmy
          sx =  sign(one,anrmx)
          sy =  sign(one,anrmy)
       else
          dg = dx_eb / abs(anrmy)
          gx = bctx - dg*anrmx
          gy = bcty - dg*anrmy
          sx =  sign(one,anrmx)
          sy =  sign(one,anrmy)
       end if
       ii = i - int(sx)
       jj = j - int(sy)
               
       w1 = amrex_blend_beta(vf)
       w2 = one-w1

       if (w1.eq.zero) then
          phig1 = zero
       else
          phig1 = (one + gx*sx + gy*sy + gx*gy*sx*sy) * phi(i,j) &
               +  (    - gx*sx         - gx*gy*sx*sy) * phi(ii,j) &
               +  (            - gy*sy - gx*gy*sx*sy) * phi(i,jj) &
               +  (                    + gx*gy*sx*sy) * phi(ii,jj) 
       end if

       if (w2.eq.zero) then
          phig2 = zero
       else
          bsxinv = one/(bctx+sx)
          bsyinv = one/(bcty+sy)
          
          c_0(0,0) = sx*sy*bsxinv*bsyinv
          c_0(-1,0) = bctx*bsxinv
          c_0(0,-1) = bcty*bsyinv
          c_0(-1,-1) = -bctx*bcty*bsxinv*bsyinv

          c_x(0,0) = sy*bsxinv*bsyinv
          c_x(-1,0) = -bsxinv
          c_x(0,-1) = sx*bcty*bsyinv
          c_x(-1,-1) = -sx*bctx*bcty*bsxinv*bsyinv

          c_y(0,0) = sx*bsxinv*bsyinv
          c_y(-1,0) = sy*bctx*bsxinv
          c_y(0,-1) = -bsyinv
          c_y(-1,-1) = -sy*bctx*bcty*bsxinv*bsyinv

          c_xy(0,0) = bsxinv*bsyinv
          c_xy(-1,0) = -sy*bsxinv
          c_xy(0,-1) = -sx*bsyinv
          c_xy(-1,-1) = (one+sx*bctx+sy*bcty)*bsxinv*bsyinv

          phig2 = (c_0( 0, 0) + gx*c_x( 0, 0) + gy*c_y( 0, 0) + gx*gy*c_xy( 0, 0)) * phib &
               +  (c_0(-1, 0) + gx*c_x(-1, 0) + gy*c_y(-1, 0) + gx*gy*c_xy(-1, 0)) * phi(ii,j) &
               +  (c_0( 0,-1) + gx*c_x( 0,-1) + gy*c_y( 0,-1) + gx*gy*c_xy( 0,-1)) * phi(i,jj) &
               +  (c_0(-1,-1) + gx*c_x(-1,-1) + gy*c_y(-1,-1) + gx*gy*c_xy(-1,-1)) * phi(ii,jj)
       end if

       phig = w1*phig1 + w2*phig2
       dphidn = (phib-phig) / dg  

  end subroutine compute_dphidn_2d

  subroutine compute_dphidn_2d_ho (dphidn, dxinv, i, j,  &
        phi , p_lo, p_hi,     &
        flag,  flo,  fhi,     &
        bct, phib, anrmx, anrmy)

      ! Cell indices 
      integer, intent(in   ) :: i, j

      ! Grid spacing
      real(amrex_real),       intent(in   ) :: dxinv(2)

      integer, intent(in   ) :: p_lo(2), p_hi(2)
      integer, intent(in   ) ::  flo(2),  fhi(2)
      integer, intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2) )

      ! Arrays
      real(amrex_real),  intent(in   ) ::                            &
           & phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2))

      real(amrex_real),  intent(in   ) :: bct(2), phib
      real(amrex_real),  intent(in   ) :: anrmx, anrmy

      real(amrex_real),        intent(  out) :: dphidn

      ! Local variable
      real(amrex_real) :: anrm
      real(amrex_real) :: xit, yit, s, s2
      real(amrex_real) :: d1, d2, ddinv
      real(amrex_real) :: u1, u2
      integer          :: ixit, iyit, is, is2, ivar

      if (abs(anrmx).ge.abs(anrmy)) then
         anrm = anrmx
         ivar = 1
      else 
         anrm = anrmy
         ivar = 2
      end if

      ! s is -1. or 1.
        s = sign(one,-anrm)
       s2 = 2*s

      ! is is -1 or 1
      is  = nint(s)
      is2 = 2*is

      d1 = (bct(ivar) - s ) * (one/anrm)
      d2 = (bct(ivar) - s2) * (one/anrm)

      if (abs(anrmx).ge.abs(anrmy)) then
         ! the equation for the line:  x = bct(1) - d*anrmx
         !                             y = bct(2) - d*anrmy

         !
         ! the line intersects (x = s) at ...
         !
         yit = bct(2) - d1*anrmy
         iyit = j   + nint(yit)
         yit  = yit - nint(yit)

         call interp1d(u1,yit,phi(i+is,iyit-1:iyit+1),flag(i+is,iyit-1:iyit+1))

         !
         ! the line intersects (x = 2*s) at ...
         !
         yit = bct(2) - d2*anrmy
         iyit = j   + nint(yit)
         yit  = yit - nint(yit)

         call interp1d(u2,yit,phi(i+is2,iyit-1:iyit+1),flag(i+is2,iyit-1:iyit+1))

      else

         xit = bct(1) - d1*anrmx
         ixit = i + nint(xit)
         xit = xit - nint(xit)

         call interp1d(u1,xit,phi(ixit-1:ixit+1,j+is),flag(ixit-1:ixit+1,j+is))

         xit = bct(1) - d2*anrmx
         ixit = i + nint(xit)
         xit = xit - nint(xit)

         call interp1d(u2,xit,phi(ixit-1:ixit+1,j+is2),flag(ixit-1:ixit+1,j+is2))

      end if

      !
      ! compute derivatives on the wall given that phi is zero on the wall.
      !
      ddinv = one/(d1*d2*(d2-d1))
      dphidn = -ddinv*( d2*d2*(u1-phib) - d1*d1*(u2-phib) )  ! note that the normal vector points toward the wall

  end subroutine compute_dphidn_2d_ho

  subroutine interp1d(phi_interp,yit,v,flag)

      use amrex_ebcellflag_module, only : is_covered_cell

      real(amrex_real), intent(  out) :: phi_interp
      real(amrex_real), intent(in   ) :: yit,v(3)
      integer         , intent(in   ) :: flag(3)

      real(amrex_real)             :: cym,cy0,cyp

      if (is_covered_cell(flag(1))) then
         phi_interp =             (one-yit)*v(2) + yit*v(3)
      else if (is_covered_cell(flag(3))) then
         phi_interp = -yit*v(1) + (one+yit)*v(2)
      else 
         ! Coefficents for quadratic interpolation
         cym = half*yit*(yit-one)
         cy0 = one-yit*yit
         cyp = half*yit*(yit+one)
         phi_interp =  cym*v(1) +       cy0*v(2) + cyp*v(3)
      end if

  end subroutine interp1d

end module amrex_mlebabeclap_2d_module
