module amrex_mlebabeclap_3d_module
  
  use amrex_error_module 
  use amrex_constants_module, only: zero, one, third, half
  use amrex_fort_module, only : amrex_real
  use amrex_ebcellflag_module, only: is_regular_cell, is_covered_cell, is_single_valued_cell, &
       get_neighbor_cells_int_single
  implicit none

  private
  public :: amrex_mlebabeclap_adotx, amrex_mlebabeclap_gsrb, amrex_mlebabeclap_normalize, & 
            amrex_eb_mg_interp, amrex_mlebabeclap_grad, amrex_mlebabeclap_flux, &
            compute_dphidn_3d, compute_dphidn_3d_ho, amrex_get_dx_eb

contains

  elemental function amrex_get_dx_eb (kappa)
    real(amrex_real), intent(in) :: kappa
    real(amrex_real) :: amrex_get_dx_eb
    amrex_get_dx_eb = max(0.3d0, (kappa*kappa-0.25d0)/(2.d0*kappa))
  end function amrex_get_dx_eb

  subroutine amrex_mlebabeclap_adotx(lo, hi, y, ylo, yhi, x, xlo, xhi, & 
       a, alo, ahi, bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, ccm, cmlo, cmhi, flag, flo, fhi, & 
       vfrc, vlo, vhi, apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi, fcx, cxlo, cxhi, &
       fcy, cylo, cyhi, fcz, czlo, czhi, ba, balo, bahi, bc, bclo, bchi, beb, elo, ehi, &
       is_dirichlet, is_ho_dirichlet, phieb, plo, phi, is_inhomog, dxinv, alpha, beta) & 
       bind(c, name='amrex_mlebabeclap_adotx') 

   integer, dimension(3), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, alo, ahi, bxlo,&
        bxhi, bylo, byhi, bzlo, bzhi, cmlo, cmhi, flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, &
        azlo, azhi, cxlo, cxhi, cylo, cyhi, czlo, czhi, balo, bahi, bclo, bchi, elo, ehi, plo, phi

   real(amrex_real), intent(in) :: dxinv(3) 
   integer         , value, intent(in) :: is_dirichlet, is_ho_dirichlet, is_inhomog
   real(amrex_real), value, intent(in) :: alpha, beta
   real(amrex_real), intent(inout) ::    y( ylo(1): yhi(1), ylo(2): yhi(2), ylo(3): yhi(3))
   real(amrex_real), intent(in   ) ::    x( xlo(1): xhi(1), xlo(2): xhi(2), xlo(3): xhi(3))
   real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
   real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
   real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
   real(amrex_real), intent(in   ) ::   bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
   integer         , intent(in   ) ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3)) 
   integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3)) 
   real(amrex_real), intent(in   ) :: vfrc( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3)) 
   real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)) 
   real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
   real(amrex_real), intent(in   ) ::  apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
   real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
   real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2) 
   real(amrex_real), intent(in   ) ::  fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2) 
   real(amrex_real), intent(in   ) ::  ba (balo(1):bahi(1),balo(2):bahi(2),balo(3):bahi(3))
   real(amrex_real), intent(in   ) ::  bc (bclo(1):bchi(1),bclo(2):bchi(2),bclo(3):bchi(3),3)
   real(amrex_real), intent(in   ) ::  beb( elo(1): ehi(1), elo(2): ehi(2), elo(3): ehi(3))
   real(amrex_real), intent(in   ) ::phieb( plo(1): phi(1), plo(2): phi(2), plo(3): phi(3))
   integer  :: i, j, k, ii, jj, kk 
   real(amrex_real) :: dhx, dhy, dhz, fxm, fxp, fym, fyp, fzm, fzp, fracx, fracy, fracz
   real(amrex_real) :: phib, feb
   real(amrex_real) :: anrmx, anrmy, anrmz, anorm, anorminv
   logical :: is_inhomogeneous

   real(amrex_real) :: dphidn

   is_inhomogeneous = is_inhomog .ne. 0

   dhx = beta*dxinv(1)*dxinv(1) 
   dhy = beta*dxinv(2)*dxinv(2)  
   dhz = beta*dxinv(3)*dxinv(3) 

   do         k = lo(3), hi(3) 
       do     j = lo(2), hi(2) 
           do i = lo(1), hi(1) 

             if(is_covered_cell(flag(i,j,k))) then 

                 y(i,j,k) = zero

              else if (is_regular_cell(flag(i,j,k))) then 

                 y(i,j,k) = alpha*a(i,j,k)*x(i,j,k) & 
                            - dhx * (bX(i+1,j,k)*(x(i+1,j,k) - x(i  ,j,k))  & 
                            &       -bX(i  ,j,k)*(x(i  ,j,k) - x(i-1,j,k))) & 
                            - dhy * (bY(i,j+1,k)*(x(i,j+1,k) - x(i,j  ,k))  &
                            &       -bY(i,j  ,k)*(x(i,j  ,k) - x(i,j-1,k))) &
                            - dhz * (bZ(i,j,k+1)*(x(i,j,k+1) - x(i,j,k  ))  & 
                            &       -bZ(i,j,k  )*(x(i,j,k  ) - x(i,j,k-1)))
                    
              else 
                fxm = bX(i,j,k)*(x(i,j,k) - x(i-1,j,k))
                if (apx(i,j,k).ne.zero.and.apx(i,j,k).ne.one) then 
                    jj = j + int(sign(one, fcx(i,j,k,1)))
                    kk = k + int(sign(one, fcx(i,j,k,2)))
                    fracy = abs(fcx(i,j,k,1))*real(ior(ccm(i-1,jj,k),ccm(i,jj,k)),amrex_real)
                    fracz = abs(fcx(i,j,k,2))*real(ior(ccm(i-1,j,kk),ccm(i,j,kk)),amrex_real)
                    fxm = (one-fracy)*(one-fracz)*fxm + &
                         & fracy*(one-fracz)*bX(i,jj,k )*(x(i,jj,k )-x(i-1,jj,k )) + & 
                         & fracz*(one-fracy)*bX(i,j ,kk)*(x(i,j ,kk)-x(i-1,j ,kk)) + &
                         & fracy*     fracz *bX(i,jj,kk)*(x(i,jj,kk)-x(i-1,jj,kk))
                endif 

                fxp = bX(i+1,j,k)*(x(i+1,j,k) - x(i,j,k))
                if (apx(i+1,j,k).ne.zero.and.apx(i+1,j,k).ne.one) then 
                    jj = j + int(sign(one,fcx(i+1,j,k,1)))
                    kk = k + int(sign(one,fcx(i+1,j,k,2)))
                    fracy = abs(fcx(i+1,j,k,1))*real(ior(ccm(i,jj,k),ccm(i+1,jj,k)),amrex_real)
                    fracz = abs(fcx(i+1,j,k,2))*real(ior(ccm(i,j,kk),ccm(i+1,j,kk)),amrex_real)
                    fxp = (one-fracy)*(one-fracz)*fxp + &
                         & fracy*(one-fracz)*bX(i+1,jj,k )*(x(i+1,jj,k )-x(i,jj,k )) + & 
                         & fracz*(one-fracy)*bX(i+1,j ,kk)*(x(i+1,j ,kk)-x(i,j ,kk)) + & 
                         & fracy*     fracz *bX(i+1,jj,kk)*(x(i+1,jj,kk)-x(i,jj,kk))
                endif 

                fym = bY(i,j,k)*(x(i,j,k) - x(i,j-1,k))
                if (apy(i,j,k).ne.zero.and.apy(i,j,k).ne.one) then 
                    ii = i + int(sign(one,fcy(i,j,k,1)))
                    kk = k + int(sign(one,fcy(i,j,k,2)))
                    fracx = abs(fcy(i,j,k,1))*real(ior(ccm(ii,j-1,k),ccm(ii,j,k)),amrex_real)
                    fracz = abs(fcy(i,j,k,2))*real(ior(ccm(i,j-1,kk),ccm(i,j,kk)),amrex_real)
                    fym = (one-fracx)*(one-fracz)*fym + &
                         & fracx*(one-fracz)*bY(ii,j,k )*(x(ii,j,k )-x(ii,j-1,k )) + & 
                         & fracz*(one-fracx)*bY(i ,j,kk)*(x(i ,j,kk)-x(i ,j-1,kk)) + &
                         & fracx*     fracz *bY(ii,j,kk)*(x(ii,j,kk)-x(ii,j-1,kk))
                endif 

                fyp = bY(i,j+1,k)*(x(i,j+1,k) - x(i,j,k))
                if (apy(i,j+1,k).ne.zero.and.apy(i,j+1,k).ne.one) then 
                    ii = i + int(sign(one,fcy(i,j+1,k,1)))
                    kk = k + int(sign(one,fcy(i,j+1,k,2)))
                    fracx = abs(fcy(i,j+1,k,1))*real(ior(ccm(ii,j,k),ccm(ii,j+1,k)),amrex_real)
                    fracz = abs(fcy(i,j+1,k,2))*real(ior(ccm(i,j,kk),ccm(i,j+1,kk)),amrex_real)
                    fyp = (one-fracx)*(one-fracz)*fyp + &
                         & fracx*(one-fracz)*bY(ii,j+1,k )*(x(ii,j+1,k )-x(ii,j,k )) + &
                         & fracz*(one-fracx)*bY(i ,j+1,kk)*(x(i ,j+1,kk)-x(i ,j,kk)) + & 
                         & fracx*     fracz *bY(ii,j+1,kk)*(x(ii,j+1,kk)-x(ii,j,kk))
                endif 

                fzm = bZ(i,j,k)*(x(i,j,k) - x(i,j,k-1))
                if (apz(i,j,k).ne.zero.and.apz(i,j,k).ne.one) then 
                    ii = i + int(sign(one,fcz(i,j,k,1)))
                    jj = j + int(sign(one,fcz(i,j,k,2)))
                    fracx = abs(fcz(i,j,k,1))*real(ior(ccm(ii,j,k-1),ccm(ii,j,k)),amrex_real)
                    fracy = abs(fcz(i,j,k,2))*real(ior(ccm(i,jj,k-1),ccm(i,jj,k)),amrex_real)
                    fzm = (one-fracx)*(one-fracy)*fzm + &
                         & fracx*(one-fracy)*bZ(ii,j ,k)*(x(ii,j ,k)-x(ii,j ,k-1)) + & 
                         & fracy*(one-fracx)*bZ(i ,jj,k)*(x(i ,jj,k)-x(i ,jj,k-1)) + &
                         & fracx*     fracy *bZ(ii,jj,k)*(x(ii,jj,k)-x(ii,jj,k-1))
                endif

                fzp = bZ(i,j,k+1)*(x(i,j,k+1) - x(i,j,k))
                if (apz(i,j,k+1).ne.zero.and.apz(i,j,k+1).ne.one) then 
                    ii = i + int(sign(one,fcz(i,j,k+1,1)))
                    jj = j + int(sign(one,fcz(i,j,k+1,2)))
                    fracx = abs(fcz(i,j,k+1,1))*real(ior(ccm(ii,j,k),ccm(ii,j,k+1)),amrex_real)
                    fracy = abs(fcz(i,j,k+1,2))*real(ior(ccm(i,jj,k),ccm(i,jj,k+1)),amrex_real)
                    fzp = (one-fracx)*(one-fracy)*fzp + & 
                         & fracx*(one-fracy)*bZ(ii,j ,k+1)*(x(ii,j ,k+1)-x(ii,j ,k)) + &
                         & fracy*(one-fracx)*bZ(i ,jj,k+1)*(x(i ,jj,k+1)-x(i ,jj,k)) + &
                         & fracx*     fracy *bZ(ii,jj,k+1)*(x(ii,jj,k+1)-x(ii,jj,k))
                endif 

                if (is_ho_dirichlet .ne. 0 .or. is_dirichlet .ne. 0)  then
                   anorm = sqrt((apx(i,j,k)-apx(i+1,j,k))**2 &
                        +       (apy(i,j,k)-apy(i,j+1,k))**2 &
                        +       (apz(i,j,k)-apz(i,j,k+1))**2)
                   anorminv = one/anorm
                   anrmx = (apx(i,j,k)-apx(i+1,j,k)) * anorminv
                   anrmy = (apy(i,j,k)-apy(i,j+1,k)) * anorminv
                   anrmz = (apz(i,j,k)-apz(i,j,k+1)) * anorminv

                   if (is_inhomogeneous) then
                      phib = phieb(i,j,k)
                   else
                      phib = zero
                   end if

                   if (is_ho_dirichlet .ne. 0) then
                      call compute_dphidn_3d_ho(dphidn, dxinv, i, j, k, &
                                                x,    xlo,  xhi, &
                                                flag, flo,  fhi, &
                                                bc(i,j,k,:),  phib, &
                                                anrmx, anrmy, anrmz)
                   else if (is_dirichlet .ne. 0) then
                      call compute_dphidn_3d(dphidn, dxinv, i, j, k, &
                                             x,    xlo,  xhi, &
                                             flag, flo,  fhi, &
                                             bc(i,j,k,:),  phib, &
                                             anrmx, anrmy, anrmz, vfrc(i,j,k))
                   end if
                   feb = dphidn * ba(i,j,k) * beb(i,j,k)
                else
                   feb = zero
                end if
                   
                y(i,j,k) = alpha*a(i,j,k)*x(i,j,k) + (one/vfrc(i,j,k))*&
                       (dhx*(apx(i,j,k)*fxm - apx(i+1,j,k)*fxp) + &
                        dhy*(apy(i,j,k)*fym - apy(i,j+1,k)*fyp) + &
                        dhz*(apz(i,j,k)*fzm - apz(i,j,k+1)*fzp) &
                        - dhx*feb)
              endif
          enddo
       enddo
    enddo
  end subroutine amrex_mlebabeclap_adotx

  subroutine amrex_mlebabeclap_gsrb(lo, hi, phi, hlo, hhi, rhs, rlo, rhi, a, alo, ahi, & 
     bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, & 
     ccm, cmlo, cmhi, &
     m0, m0lo, m0hi, m2, m2lo, m2hi, m4, m4lo, m4hi, & 
     m1, m1lo, m1hi, m3, m3lo, m3hi, m5, m5lo, m5hi, & 
     f0, f0lo, f0hi, f2, f2lo, f2hi, f4, f4lo, f4hi, & 
     f1, f1lo, f1hi, f3, f3lo, f3hi, f5, f5lo, f5hi, & 
     flag, flo, fhi, vfrc, vlo, vhi, & 
     apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
     fcz, czlo, czhi, ba, balo, bahi, bc, bclo, bchi, beb, elo, ehi, &
     is_dirichlet, is_ho_dirichlet, &
     dxinv, alpha, beta, redblack) & 
     bind(c,name='amrex_mlebabeclap_gsrb') 

    integer, dimension(3), intent(in) :: lo, hi, hlo, hhi, rlo, rhi, alo, ahi, bxlo, bxhi, bylo, byhi, &
         bzlo, bzhi, cmlo, cmhi, &
         m0lo, m0hi, m1lo, m1hi, m2lo, m2hi, m3lo, m3hi, m4lo, m4hi, m5lo, m5hi,  &
         f0lo, f0hi, f1lo, f1hi, f2lo, f2hi, f3lo, f3hi, f4lo, f4hi ,f5lo, f5hi, flo, fhi, vlo, vhi,&
         axlo, axhi, aylo, ayhi, azlo, azhi, cxlo, cxhi, cylo, cyhi, czlo, czhi, &
         balo, bahi, bclo, bchi, elo, ehi
    real(amrex_real), intent(in) :: dxinv(3)
    integer         , value, intent(in) :: is_dirichlet, is_ho_dirichlet
    real(amrex_real), value, intent(in) :: alpha, beta
    integer         , value, intent(in) :: redblack
    real(amrex_real), intent(inout) ::  phi( hlo(1): hhi(1), hlo(2): hhi(2), hlo(3): hhi(3)  )
    real(amrex_real), intent(in   ) ::  rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3)  )
    real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3)  )
    real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3)  )
    real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3)  )
    real(amrex_real), intent(in   ) ::   bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3)  )
    integer         , intent(in   ) ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3)  ) 
    integer         , intent(in   ) ::   m0(m0lo(1):m0hi(1),m0lo(2):m0hi(2),m0lo(3):m0hi(3)  )
    integer         , intent(in   ) ::   m1(m1lo(1):m1hi(1),m1lo(2):m1hi(2),m1lo(3):m1hi(3)  )
    integer         , intent(in   ) ::   m2(m2lo(1):m2hi(1),m2lo(2):m2hi(2),m2lo(3):m2hi(3)  )
    integer         , intent(in   ) ::   m3(m3lo(1):m3hi(1),m3lo(2):m3hi(2),m3lo(3):m3hi(3)  )
    integer         , intent(in   ) ::   m4(m4lo(1):m4hi(1),m4lo(2):m4hi(2),m4lo(3):m4hi(3)  )
    integer         , intent(in   ) ::   m5(m5lo(1):m5hi(1),m5lo(2):m5hi(2),m5lo(3):m5hi(3)  )
    real(amrex_real), intent(in   ) ::   f0(f0lo(1):f0hi(1),f0lo(2):f0hi(2),f0lo(3):f0hi(3)  )
    real(amrex_real), intent(in   ) ::   f1(f1lo(1):f1hi(1),f1lo(2):f1hi(2),f1lo(3):f1hi(3)  )
    real(amrex_real), intent(in   ) ::   f2(f2lo(1):f2hi(1),f2lo(2):f2hi(2),f2lo(3):f2hi(3)  )
    real(amrex_real), intent(in   ) ::   f3(f3lo(1):f3hi(1),f3lo(2):f3hi(2),f3lo(3):f3hi(3)  )
    real(amrex_real), intent(in   ) ::   f4(f4lo(1):f4hi(1),f4lo(2):f4hi(2),f4lo(3):f4hi(3)  )
    real(amrex_real), intent(in   ) ::   f5(f5lo(1):f5hi(1),f5lo(2):f5hi(2),f5lo(3):f5hi(3)  )
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3)  )
    real(amrex_real), intent(in   ) :: vfrc( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3)  )
    real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)  )
    real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)  )
    real(amrex_real), intent(in   ) ::  apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3)  )
    real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
    real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2)
    real(amrex_real), intent(in   ) ::  fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2)
    real(amrex_real), intent(in   ) ::  ba (balo(1):bahi(1),balo(2):bahi(2),balo(3):bahi(3))
    real(amrex_real), intent(in   ) ::  bc (bclo(1):bchi(1),bclo(2):bchi(2),bclo(3):bchi(3),3)
    real(amrex_real), intent(in   ) ::  beb( elo(1): ehi(1), elo(2): ehi(2), elo(3): ehi(3))

    integer :: i, j, k, ioff, ii, jj, kk
    real(amrex_real) :: cf0, cf1, cf2, cf3, cf4, cf5,  delta, gamma, rho, res, vfrcinv
    real(amrex_real) :: dhx, dhy, dhz, fxm, fxp, fym, fyp, fzm, fzp, fracx, fracy, fracz
    real(amrex_real) :: sxm, sxp, sym, syp, szm, szp, oxm, oxp, oym, oyp, ozm, ozp
    real(amrex_real) :: feb, phig, gx, gy, gz, dg, dx_eb, gxy, gxz, gyz, gxyz
    real(amrex_real) :: feb_gamma, phig_gamma
    real(amrex_real) :: anrmx, anrmy, anrmz, anorm, anorminv, sx, sy, sz
    real(amrex_real) :: bctx, bcty, bctz
    real(amrex_real), parameter :: omega = 1.15_amrex_real ! over-relaxation

    real(amrex_real) :: dphidn, phib

    dhx = beta*dxinv(1)*dxinv(1) 
    dhy = beta*dxinv(2)*dxinv(2) 
    dhz = beta*dxinv(3)*dxinv(3) 

    do       k = lo(3), hi(3) 
       do    j = lo(2), hi(2)
          ioff = mod(lo(1)+k+j+redblack,2)
          do i = lo(1)+ioff, hi(1), 2
             if (is_covered_cell(flag(i,j,k))) then 
                phi(i,j,k) = zero 
             else 
                cf0 = merge(f0(lo(1),j,k), 0.0D0, & 
                      (i .eq. lo(1)) .and. (m0(lo(1)-1,j,k).gt.0))
                cf1 = merge(f1(i,lo(2),k), 0.0D0, & 
                      (j .eq. lo(2)) .and. (m1(i,lo(2)-1,k).gt.0)) 
                cf2 = merge(f2(i,j,lo(3)), 0.0D0, & 
                      (k .eq. lo(3)) .and. (m2(i,j,lo(3)-1).gt.0))
                cf3 = merge(f3(hi(1),j,k), 0.0D0, & 
                      (i .eq. hi(1)) .and. (m3(hi(1)+1,j,k).gt.0))
                cf4 = merge(f4(i,hi(2),k), 0.0D0, & 
                      (j .eq. hi(2)) .and. (m4(i,hi(2)+1,k).gt.0))
                cf5 = merge(f5(i,j,hi(3)), 0.0D0, & 
                      (k .eq. hi(3)) .and. (m5(i,j,hi(3)+1).gt.0)) 

                if (is_regular_cell(flag(i,j,k))) then 
                  
                   gamma = alpha*a(i,j,k) & 
                         + dhx*(bX(i+1,j,k) + bX(i,j,k)) & 
                         + dhy*(bY(i,j+1,k) + bY(i,j,k)) &
                         + dhz*(bZ(i,j,k+1) + bZ(i,j,k)) 

                   rho   = dhx*(bX(i+1,j,k)*phi(i+1,j,k) + bX(i,j,k)*phi(i-1,j,k)) & 
                         + dhy*(bY(i,j+1,k)*phi(i,j+1,k) + bY(i,j,k)*phi(i,j-1,k)) & 
                         + dhz*(bZ(i,j,k+1)*phi(i,j,k+1) + bZ(i,j,k)*phi(i,j,k-1))

                   delta = dhx*(bX(i,j,k)*cf0 + bX(i+1,j,k)*cf3) & 
                        +  dhy*(bY(i,j,k)*cf1 + bY(i,j+1,k)*cf4) & 
                        +  dhz*(bZ(i,j,k)*cf2 + bZ(i,j,k+1)*cf5)

                else 

                   fxm = -bX(i,j,k)*phi(i-1,j,k)
                   oxm = -bX(i,j,k)*cf0
                   sxm =  bX(i,j,k) 
                   if(apx(i,j,k).ne.zero .and. apx(i,j,k).ne.one) then 
                      jj = j + int(sign(one, fcx(i,j,k,1)))
                      kk = k + int(sign(one, fcx(i,j,k,2)))
                      fracy = abs(fcx(i,j,k,1))*real(ior(ccm(i-1,jj,k),ccm(i,jj,k)),amrex_real)
                      fracz = abs(fcx(i,j,k,2))*real(ior(ccm(i-1,j,kk),ccm(i,j,kk)),amrex_real)
                      fxm = (one-fracy)*(one-fracz)*fxm &
                           +     fracy *(one-fracz)*bX(i,jj,k )*(phi(i,jj,k )-phi(i-1,jj,k )) & 
                           +(one-fracy)*     fracz *bX(i,j ,kk)*(phi(i,j ,kk)-phi(i-1,j ,kk)) &
                           +     fracy *     fracz *bX(i,jj,kk)*(phi(i,jj,kk)-phi(i-1,jj,kk))
                      ! oxm = (one-fracy)*(one-fracz)*oxm
                      oxm = zero
                      sxm = (one-fracy)*(one-fracz)*sxm
                   end if
                   
                   fxp =  bX(i+1,j,k)*phi(i+1,j,k)
                   oxp =  bX(i+1,j,k)*cf3
                   sxp = -bX(i+1,j,k)
                   if(apx(i+1,j,k).ne.zero.and.apx(i+1,j,k).ne.one) then 
                      jj = j + int(sign(one, fcx(i+1,j,k,1)))
                      kk = k + int(sign(one, fcx(i+1,j,k,2)))
                      fracy = abs(fcx(i+1,j,k,1))*real(ior(ccm(i,jj,k),ccm(i+1,jj,k)),amrex_real)
                      fracz = abs(fcx(i+1,j,k,2))*real(ior(ccm(i,j,kk),ccm(i+1,j,kk)),amrex_real)
                      fxp = (one-fracy)*(one-fracz)*fxp &
                           +     fracy *(one-fracz)*bX(i+1,jj,k )*(phi(i+1,jj,k )-phi(i,jj,k )) & 
                           +(one-fracy)*     fracz *bX(i+1,j ,kk)*(phi(i+1,j ,kk)-phi(i,j ,kk)) & 
                           +     fracy *     fracz *bX(i+1,jj,kk)*(phi(i+1,jj,kk)-phi(i,jj,kk))
                      ! oxp = (one-fracy)*(one-fracz)*oxp
                      oxp = zero
                      sxp = (one-fracy)*(one-fracz)*sxp
                   end if
                   
                   fym = -bY(i,j,k)*phi(i,j-1,k)
                   oym = -bY(i,j,k)*cf1
                   sym =  bY(i,j,k)
                   if(apy(i,j,k).ne.zero.and.apy(i,j,k).ne.one) then 
                      ii = i + int(sign(one,fcy(i,j,k,1)))
                      kk = k + int(sign(one,fcy(i,j,k,2)))
                      fracx = abs(fcy(i,j,k,1))*real(ior(ccm(ii,j-1,k),ccm(ii,j,k)),amrex_real)
                      fracz = abs(fcy(i,j,k,2))*real(ior(ccm(i,j-1,kk),ccm(i,j,kk)),amrex_real)
                      fym = (one-fracx)*(one-fracz)*fym &
                           +     fracx *(one-fracz)*bY(ii,j,k )*(phi(ii,j,k )-phi(ii,j-1,k )) & 
                           +(one-fracx)*     fracz *bY(i ,j,kk)*(phi(i ,j,kk)-phi(i ,j-1,kk)) &
                           +     fracx *     fracz *bY(ii,j,kk)*(phi(ii,j,kk)-phi(ii,j-1,kk))
                      ! oym = (one-fracx)*(one-fracz)*oym
                      oym = zero
                      sym = (one-fracx)*(one-fracz)*sym
                   endif
 
                   fyp =  bY(i,j+1,k)*phi(i,j+1,k)
                   oyp =  bY(i,j+1,k)*cf4
                   syp = -bY(i,j+1,k)
                   if(apy(i,j+1,k).ne.zero.and.apy(i,j+1,k).ne.one) then 
                      ii = i + int(sign(one,fcy(i,j+1,k,1)))
                      kk = k + int(sign(one,fcy(i,j+1,k,2)))
                      fracx = abs(fcy(i,j+1,k,1))*real(ior(ccm(ii,j,k),ccm(ii,j+1,k)),amrex_real)
                      fracz = abs(fcy(i,j+1,k,2))*real(ior(ccm(i,j,kk),ccm(i,j+1,kk)),amrex_real)
                      fyp = (one-fracx)*(one-fracz)*fyp &
                           +     fracx *(one-fracz)*bY(ii,j+1,k )*(phi(ii,j+1,k )-phi(ii,j,k )) &
                           +(one-fracx)*     fracz *bY(i ,j+1,kk)*(phi(i ,j+1,kk)-phi(i ,j,kk)) & 
                           +     fracx *     fracz *bY(ii,j+1,kk)*(phi(ii,j+1,kk)-phi(ii,j,kk))
                      ! oyp = (one-fracx)*(one-fracz)*oyp
                      oyp = zero
                      syp = (one-fracx)*(one-fracz)*syp
                   end if
 
                   fzm = -bZ(i,j,k)*phi(i,j,k-1)
                   ozm = -bZ(i,j,k)*cf2
                   szm =  bZ(i,j,k)
                   if(apz(i,j,k).ne.zero.and.apz(i,j,k).ne.one) then 
                      ii = i + int(sign(one,fcz(i,j,k,1)))
                      jj = j + int(sign(one,fcz(i,j,k,2)))
                      fracx = abs(fcz(i,j,k,1))*real(ior(ccm(ii,j,k-1),ccm(ii,j,k)),amrex_real)
                      fracy = abs(fcz(i,j,k,2))*real(ior(ccm(i,jj,k-1),ccm(i,jj,k)),amrex_real)
                      fzm = (one-fracx)*(one-fracy)*fzm &
                           +     fracx *(one-fracy)*bZ(ii,j ,k)*(phi(ii,j ,k)-phi(ii,j ,k-1)) & 
                           +(one-fracx)*     fracy *bZ(i ,jj,k)*(phi(i ,jj,k)-phi(i ,jj,k-1)) &
                           +     fracx *     fracy *bZ(ii,jj,k)*(phi(ii,jj,k)-phi(ii,jj,k-1))
                      ! ozm = (one-fracx)*(one-fracy)*ozm
                      ozm = zero
                      szm = (one-fracx)*(one-fracy)*szm
                    endif
               
                    fzp =  bZ(i,j,k+1)*phi(i,j,k+1)
                    ozp =  bZ(i,j,k+1)*cf5
                    szp = -bZ(i,j,k+1)
                    if(apz(i,j,k+1).ne.zero.and.apz(i,j,k+1).ne.one) then 
                       ii = i + int(sign(one,fcz(i,j,k+1,1)))
                       jj = j + int(sign(one,fcz(i,j,k+1,2)))
                       fracx = abs(fcz(i,j,k+1,1))*real(ior(ccm(ii,j,k),ccm(ii,j,k+1)),amrex_real)
                       fracy = abs(fcz(i,j,k+1,2))*real(ior(ccm(i,jj,k),ccm(i,jj,k+1)),amrex_real)
                       fzp = (one-fracx)*(one-fracy)*fzp & 
                            +     fracx *(one-fracy)*bZ(ii,j ,k+1)*(phi(ii,j ,k+1)-phi(ii,j ,k)) &
                            +(one-fracx)*     fracy *bZ(i ,jj,k+1)*(phi(i ,jj,k+1)-phi(i ,jj,k)) &
                            +     fracx *     fracy *bZ(ii,jj,k+1)*(phi(ii,jj,k+1)-phi(ii,jj,k))
                       ! ozp = (one-fracx)*(one-fracy)*ozp
                       ozp = zero
                       szp = (one-fracx)*(one-fracy)*szp
                    end if 
 
                    vfrcinv = one/vfrc(i,j,k)
                    gamma = alpha*a(i,j,k) + vfrcinv * & 
                            (dhx*(apx(i,j,k)*sxm-apx(i+1,j,k)*sxp) + &
                             dhy*(apy(i,j,k)*sym-apy(i,j+1,k)*syp) + &
                             dhz*(apz(i,j,k)*szm-apz(i,j,k+1)*szp))

                    rho = -vfrcinv * & 
                           (dhx*(apx(i,j,k)*fxm-apx(i+1,j,k)*fxp) + &
                            dhy*(apy(i,j,k)*fym-apy(i,j+1,k)*fyp) + &
                            dhz*(apz(i,j,k)*fzm-apz(i,j,k+1)*fzp))

                    delta = -vfrcinv * & 
                         (dhx*(apx(i,j,k)*oxm-apx(i+1,j,k)*oxp) + &
                          dhy*(apy(i,j,k)*oym-apy(i,j+1,k)*oyp) + &
                          dhz*(apz(i,j,k)*ozm-apz(i,j,k+1)*ozp))

                    if (is_ho_dirichlet .ne. 0 .or. is_dirichlet .eq. 0) then

                       anorm = sqrt((apx(i,j,k)-apx(i+1,j,k))**2 &
                            +       (apy(i,j,k)-apy(i,j+1,k))**2 &
                            +       (apz(i,j,k)-apz(i,j,k+1))**2)
                       anorminv = one/anorm
                       anrmx = (apx(i,j,k)-apx(i+1,j,k)) * anorminv
                       anrmy = (apy(i,j,k)-apy(i,j+1,k)) * anorminv
                       anrmz = (apz(i,j,k)-apz(i,j,k+1)) * anorminv

                       ! In gsrb we are always in residual-correction form so phib = 0
                       phib = zero 

                    end if

                    if (is_ho_dirichlet .ne. 0) then
                       call compute_dphidn_3d_ho(dphidn, dxinv, i, j, k, &
                                                 phi,  hlo,  hhi, &
                                                 flag, flo,  fhi, &
                                                 bc(i,j,k,:), phib,     &
                                                 anrmx, anrmy, anrmz)
   
                       ! We should modify these but haven't done it yet
                       ! feb_gamma = -phig_gamma * (ba(i,j,k)*beb(i,j,k)/dg)
                       ! gamma = gamma + vfrcinv*(-dhx)*feb_gamma

                       feb = dphidn * ba(i,j,k) * beb(i,j,k)
                       rho = rho - vfrcinv*(-dhx)*feb

                    else if (is_dirichlet .ne. 0) then

                       anorm = sqrt((apx(i,j,k)-apx(i+1,j,k))**2 &
                            +       (apy(i,j,k)-apy(i,j+1,k))**2 &
                            +       (apz(i,j,k)-apz(i,j,k+1))**2)
                       anorminv = one/anorm
                       anrmx = (apx(i,j,k)-apx(i+1,j,k)) * anorminv
                       anrmy = (apy(i,j,k)-apy(i,j+1,k)) * anorminv
                       anrmz = (apz(i,j,k)-apz(i,j,k+1)) * anorminv
                       bctx = bc(i,j,k,1)
                       bcty = bc(i,j,k,2)
                       bctz = bc(i,j,k,3)

                       dx_eb = amrex_get_dx_eb(vfrc(i,j,k))
                       dg = dx_eb / max(abs(anrmx),abs(anrmy),abs(anrmz))

                       gx = bctx - dg*anrmx
                       gy = bcty - dg*anrmy
                       gz = bctz - dg*anrmz
                       sx =  sign(one,anrmx)
                       sy =  sign(one,anrmy)
                       sz =  sign(one,anrmz)
                       ii = i - int(sx)
                       jj = j - int(sy)
                       kk = k - int(sz)

                       gx = sx*gx
                       gy = sy*gy
                       gz = sz*gz
                       gxy = gx*gy
                       gxz = gx*gz
                       gyz = gy*gz
                       gxyz = gx*gy*gz
                       phig_gamma = (one+gx+gy+gz+gxy+gxz+gyz+gxyz)
                       phig = (-gz - gxz - gyz - gxyz) * phi(i,j,kk) &
                            + (-gy - gxy - gyz - gxyz) * phi(i,jj,k) &
                            + (gyz + gxyz) * phi(i,jj,kk) &
                            + (-gx - gxy - gxz - gxyz) * phi(ii,j,k) &
                            + (gxz + gxyz) * phi(ii,j,kk) &
                            + (gxy + gxyz) * phi(ii,jj,k) &
                            + (-gxyz) * phi(ii,jj,kk)

                       dphidn    = (    -phig)/dg
                       feb_gamma = -phig_gamma/dg * ba(i,j,k) * beb(i,j,k)
                    
                       gamma = gamma + vfrcinv*(-dhx)*feb_gamma
                       feb = dphidn * ba(i,j,k) * beb(i,j,k)
                       rho = rho - vfrcinv*(-dhx)*feb
                       
                    end if

                 end if

                 res = rhs(i,j,k) - (gamma*phi(i,j,k) - rho)
                 phi(i,j,k) = phi(i,j,k) + omega*res/(gamma-delta)

              endif
          end do 
       end do 
    enddo        

  end subroutine amrex_mlebabeclap_gsrb

  subroutine amrex_mlebabeclap_normalize (lo, hi, x, xlo, xhi, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, ccm, cmlo, cmhi, flag, flo, fhi, vfrc, vlo, vhi, &
       apx, axlo, axhi, apy, aylo, ayhi,apz, azlo, azhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
       fcz, czlo, czhi, ba, balo, bahi, bc, bclo, bchi, beb, elo, ehi, &
       is_dirichlet, is_ho_dirichlet, dxinv, alpha, beta) &
       bind(c,name='amrex_mlebabeclap_normalize')

    integer, dimension(3), intent(in) :: lo, hi, xlo, xhi, alo, ahi, bxlo, bxhi, bylo, byhi, bzlo, bzhi, &
         cmlo, cmhi, flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi,azlo, azhi, &
         cxlo, cxhi, cylo, cyhi, czlo, czhi, balo, bahi, bclo, bchi, elo, ehi
    integer         , value, intent(in) :: is_dirichlet, is_ho_dirichlet
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), value, intent(in) :: alpha, beta
    real(amrex_real), intent(inout) ::    x( xlo(1): xhi(1), xlo(2): xhi(2), xlo(3): xhi(3)  )
    real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3)  )
    real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bzhi(3)  )
    real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3)  )
    real(amrex_real), intent(in   ) ::   bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3)  )
    integer         , intent(in   ) ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3)  ) 
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3)  )
    real(amrex_real), intent(in   ) :: vfrc( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3)  )
    real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)  )
    real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)  )
    real(amrex_real), intent(in   ) ::  apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3)  )
    real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
    real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2)
    real(amrex_real), intent(in   ) ::  fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2)
    real(amrex_real), intent(in   ) ::  ba (balo(1):bahi(1),balo(2):bahi(2),balo(3):bahi(3))
    real(amrex_real), intent(in   ) ::  bc (bclo(1):bchi(1),bclo(2):bchi(2),bclo(3):bchi(3),3)
    real(amrex_real), intent(in   ) ::  beb( elo(1): ehi(1), elo(2): ehi(2), elo(3): ehi(3))

    integer :: i, j, k, ii, jj, kk
    real(amrex_real) :: dhx, dhy, dhz, sxm, sxp, sym, syp, szm, szp, gamma, fracx, fracy, fracz, vfrcinv
    real(amrex_real) :: gx, gy, gz, dg, dx_eb, gxy, gxz, gyz, gxyz
    real(amrex_real) :: feb_gamma, phig_gamma
    real(amrex_real) :: anrmx, anrmy, anrmz, anorm, anorminv, sx, sy, sz
    real(amrex_real) :: bctx, bcty, bctz

    dhx = beta*dxinv(1)*dxinv(1)
    dhy = beta*dxinv(2)*dxinv(2)
    dhz = beta*dxinv(3)*dxinv(3) 

    do      k = lo(3), hi(3)
      do    j = lo(2), hi(2)
        do  i = lo(1), hi(1)
          if (is_regular_cell(flag(i,j,k))) then
             x(i,j,k) = x(i,j,k) / (alpha*a(i,j,k) + dhx*(bX(i,j,k)+bX(i+1,j,k)) &
                  &                                + dhy*(bY(i,j,k)+bY(i,j+1,k)) &
                  &                                + dhz*(bZ(i,j,k)+bZ(i,j,k+1)))
          else if (is_single_valued_cell(flag(i,j,k))) then

             sxm =  bX(i,j,k) 
             if(apx(i,j,k).ne.zero .and. apx(i,j,k).ne.one) then 
                jj = j + int(sign(one, fcx(i,j,k,1)))
                kk = k + int(sign(one, fcx(i,j,k,2)))
                fracy = abs(fcx(i,j,k,1))*real(ior(ccm(i-1,jj,k),ccm(i,jj,k)),amrex_real)
                fracz = abs(fcx(i,j,k,2))*real(ior(ccm(i-1,j,kk),ccm(i,j,kk)),amrex_real)
                sxm = (one-fracy)*(one-fracz)*sxm
             end if
                   
             sxp = -bX(i+1,j,k)
             if(apx(i+1,j,k).ne.zero.and.apx(i+1,j,k).ne.one) then 
                jj = j + int(sign(one, fcx(i+1,j,k,1)))
                kk = k + int(sign(one, fcx(i+1,j,k,2)))
                fracy = abs(fcx(i+1,j,k,1))*real(ior(ccm(i,jj,k),ccm(i+1,jj,k)),amrex_real)
                fracz = abs(fcx(i+1,j,k,2))*real(ior(ccm(i,j,kk),ccm(i+1,j,kk)),amrex_real)
                sxp = (one-fracy)*(one-fracz)*sxp
             end if
                   
             sym =  bY(i,j,k)
             if(apy(i,j,k).ne.zero.and.apy(i,j,k).ne.one) then 
                ii = i + int(sign(one,fcy(i,j,k,1)))
                kk = k + int(sign(one,fcy(i,j,k,2)))
                fracx = abs(fcy(i,j,k,1))*real(ior(ccm(ii,j-1,k),ccm(ii,j,k)),amrex_real)
                fracz = abs(fcy(i,j,k,2))*real(ior(ccm(i,j-1,kk),ccm(i,j,kk)),amrex_real)
                sym = (one-fracx)*(one-fracz)*sym
             endif
 
             syp = -bY(i,j+1,k)
             if(apy(i,j+1,k).ne.zero.and.apy(i,j+1,k).ne.one) then 
                ii = i + int(sign(one,fcy(i,j+1,k,1)))
                kk = k + int(sign(one,fcy(i,j+1,k,2)))
                fracx = abs(fcy(i,j+1,k,1))*real(ior(ccm(ii,j,k),ccm(ii,j+1,k)),amrex_real)
                fracz = abs(fcy(i,j+1,k,2))*real(ior(ccm(i,j,kk),ccm(i,j+1,kk)),amrex_real)
                syp = (one-fracx)*(one-fracz)*syp
             end if
             
             szm =  bZ(i,j,k)
             if(apz(i,j,k).ne.zero.and.apz(i,j,k).ne.one) then 
                ii = i + int(sign(one,fcz(i,j,k,1)))
                jj = j + int(sign(one,fcz(i,j,k,2)))
                fracx = abs(fcz(i,j,k,1))*real(ior(ccm(ii,j,k-1),ccm(ii,j,k)),amrex_real)
                fracy = abs(fcz(i,j,k,2))*real(ior(ccm(i,jj,k-1),ccm(i,jj,k)),amrex_real)
                szm = (one-fracx)*(one-fracy)*szm
             endif
               
             szp = -bZ(i,j,k+1)
             if(apz(i,j,k+1).ne.zero.and.apz(i,j,k+1).ne.one) then 
                ii = i + int(sign(one,fcz(i,j,k+1,1)))
                jj = j + int(sign(one,fcz(i,j,k+1,2)))
                fracx = abs(fcz(i,j,k+1,1))*real(ior(ccm(ii,j,k),ccm(ii,j,k+1)),amrex_real)
                fracy = abs(fcz(i,j,k+1,2))*real(ior(ccm(i,jj,k),ccm(i,jj,k+1)),amrex_real)
                szp = (one-fracx)*(one-fracy)*szp
             end if

             vfrcinv = one/vfrc(i,j,k)
             gamma = alpha*a(i,j,k) + vfrcinv * &
                  (dhx*(apx(i,j,k)*sxm-apx(i+1,j,k)*sxp) + &
                   dhy*(apy(i,j,k)*sym-apy(i,j+1,k)*syp) + & 
                   dhz*(apz(i,j,k)*szm-apz(i,j,k+1)*szp))

             if (is_dirichlet .ne. 0 .or. is_ho_dirichlet .ne. 0) then
                anorm = sqrt((apx(i,j,k)-apx(i+1,j,k))**2 &
                     +       (apy(i,j,k)-apy(i,j+1,k))**2 &
                     +       (apz(i,j,k)-apz(i,j,k+1))**2)
                anorminv = one/anorm
                anrmx = (apx(i,j,k)-apx(i+1,j,k)) * anorminv
                anrmy = (apy(i,j,k)-apy(i,j+1,k)) * anorminv
                anrmz = (apz(i,j,k)-apz(i,j,k+1)) * anorminv
                bctx = bc(i,j,k,1)
                bcty = bc(i,j,k,2)
                bctz = bc(i,j,k,3)

                dx_eb = amrex_get_dx_eb(vfrc(i,j,k))
                dg = dx_eb / max(abs(anrmx),abs(anrmy),abs(anrmz))

                gx = bctx - dg*anrmx
                gy = bcty - dg*anrmy
                gz = bctz - dg*anrmz
                sx =  sign(one,anrmx)
                sy =  sign(one,anrmy)
                sz =  sign(one,anrmz)
                ii = i - int(sx)
                jj = j - int(sy)
                kk = k - int(sz)
                
                gx = sx*gx
                gy = sy*gy
                gz = sz*gz
                gxy = gx*gy
                gxz = gx*gz
                gyz = gy*gz
                gxyz = gx*gy*gz
                phig_gamma = (one+gx+gy+gz+gxy+gxz+gyz+gxyz)
                feb_gamma = -phig_gamma * (ba(i,j,k)*beb(i,j,k)/dg)
                    
                gamma = gamma + vfrcinv*(-dhx)*feb_gamma
             end if

             x(i,j,k) = x(i,j,k) / gamma
          end if
        end do
      end do
    end do
  end subroutine amrex_mlebabeclap_normalize


  subroutine amrex_eb_mg_interp (lo, hi, fine, flo, fhi, crse, clo, chi, flag, glo, ghi, ncomp) &
       bind(c,name='amrex_eb_mg_interp')
    integer, dimension(3), intent(in) :: lo, hi, flo, fhi, clo, chi, glo, ghi
    integer, intent(in) :: ncomp
    real(amrex_real), intent(inout) :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),ncomp)
    real(amrex_real), intent(in   ) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),ncomp)
    integer         , intent(in   ) :: flag(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3))

    integer :: i,j,ii,jj,k,kk,n

    do n = 1, ncomp
       do k = lo(3), hi(3) 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               kk = 2*k
               ii = 2*i
               jj = 2*j
               if (.not.is_covered_cell(flag(ii,jj,kk))) then
                  fine(ii,jj,kk,n) = fine(ii,jj,kk,n) + crse(i,j,k,n)
               end if

               ii = 2*i+1
               if (.not.is_covered_cell(flag(ii,jj,kk))) then
                  fine(ii,jj,kk,n) = fine(ii,jj,kk,n) + crse(i,j,k,n)
               end if

               ii = 2*i
               jj = 2*j+1
               if (.not.is_covered_cell(flag(ii,jj,kk))) then
                  fine(ii,jj,kk,n) = fine(ii,jj,kk,n) + crse(i,j,k,n)
               end if

               ii = 2*i+1
               if (.not.is_covered_cell(flag(ii,jj,kk))) then
                  fine(ii,jj,kk,n) = fine(ii,jj,kk,n) + crse(i,j,k,n)
               end if

               kk = 2*k+1
               ii = 2*i
               jj = 2*j
               if (.not.is_covered_cell(flag(ii,jj,kk))) then
                  fine(ii,jj,kk,n) = fine(ii,jj,kk,n) + crse(i,j,k,n)
               end if

               ii = 2*i+1
               if (.not.is_covered_cell(flag(ii,jj,kk))) then
                  fine(ii,jj,kk,n) = fine(ii,jj,kk,n) + crse(i,j,k,n)
               end if

               ii = 2*i
               jj = 2*j+1
               if (.not.is_covered_cell(flag(ii,jj,kk))) then
                  fine(ii,jj,kk,n) = fine(ii,jj,kk,n) + crse(i,j,k,n)
               end if

               ii = 2*i+1
               if (.not.is_covered_cell(flag(ii,jj,kk))) then
                  fine(ii,jj,kk,n) = fine(ii,jj,kk,n) + crse(i,j,k,n)
               end if

            end do
         end do
      end do
   end do 
  end subroutine amrex_eb_mg_interp

  subroutine amrex_mlebabeclap_flux(lo, hi, fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, &
                                    apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi,      &
                                    fcx, cxlo, cxhi, fcy, cylo, cyhi, fcz, czlo, czhi,      &
                                    x, slo, shi, bx, bxlo, bxhi, by, bylo, byhi,          &
                                    bz, bzlo, bzhi, ccm, cmlo, cmhi, &
                                    flag, flo, fhi, dxinv, beta, face_only) &
                                    bind(c, name='amrex_mlebabeclap_flux')

    integer, dimension(3), intent(in) :: lo, hi, slo, shi, bxlo, bxhi, bylo, byhi, bzlo, bzhi, &
         flo, fhi, axlo, axhi, aylo, ayhi, azlo, azhi,         &
         fxlo,fxhi,fylo, fyhi, fzlo, fzhi,  cxlo, cxhi, cylo,  &
         cyhi, czlo, czhi, cmlo, cmhi
    real(amrex_real), intent(in) :: dxinv(3) 
    real(amrex_real), value, intent(in) ::  beta
    real(amrex_real), intent(inout) ::   fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
    real(amrex_real), intent(inout) ::   fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
    real(amrex_real), intent(inout) ::   fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
    real(amrex_real), intent(in   ) ::   x( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
    real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
    real(amrex_real), intent(in   ) ::   bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
    integer         , intent(in   ) ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3)) 
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3)) 
    integer, value  , intent(in   ) :: face_only 
    real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)) 
    real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(amrex_real), intent(in   ) ::  apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
    real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2) 
    real(amrex_real), intent(in   ) ::  fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2) 
    integer  :: i, j, k, ii, jj, kk, istride, jstride, kstride
    real(amrex_real) :: dhx, dhy, dhz, fxm, fym, fzm, fracx, fracy, fracz
    
    dhx = beta*dxinv(1)
    dhy = beta*dxinv(2)
    dhz = beta*dxinv(3)
    
    if (face_only .eq. 1) then
       istride = hi(1)+1-lo(1)
       jstride = hi(2)+1-lo(2)
       kstride = hi(3)+1-lo(3)
    else
       istride = 1
       jstride = 1
       kstride = 1
    end if
    
    do       k = lo(3), hi(3) 
       do    j = lo(2), hi(2) 
          do i = lo(1), hi(1)+1, istride
             if(apx(i,j,k) .eq. zero) then
                fx(i,j,k) = zero
             else if (is_regular_cell(flag(i,j,k)) .or. apx(i,j,k).eq.one) then
                fx(i,j,k) = - dhx*bX(i,j,k)*(x(i,j,k) - x(i-1,j,k)) 
             else 
                fxm = bX(i,j,k)*(x(i,j,k) - x(i-1,j,k))
                jj = j + int(sign(one, fcx(i,j,k,1)))
                kk = k + int(sign(one, fcx(i,j,k,2)))
                fracy = abs(fcx(i,j,k,1))*real(ior(ccm(i-1,jj,k),ccm(i,jj,k)),amrex_real)
                fracz = abs(fcx(i,j,k,2))*real(ior(ccm(i-1,j,kk),ccm(i,j,kk)),amrex_real)
                fxm = (one-fracy)*(one-fracz)*fxm + &
                     & fracy*(one-fracz)*bX(i,jj,k )*(x(i,jj,k )-x(i-1,jj,k )) + & 
                     & fracz*(one-fracy)*bX(i,j ,kk)*(x(i,j ,kk)-x(i-1,j ,kk)) + &
                     & fracy*     fracz *bX(i,jj,kk)*(x(i,jj,kk)-x(i-1,jj,kk))
                fx(i,j,k) = -dhx*fxm
             endif
          enddo
       enddo
    enddo
    
    do       k = lo(3), hi(3) 
       do    j = lo(2), hi(2)+1, jstride 
          do i = lo(1), hi(1) 
             if(apy(i,j,k) .eq. zero) then
                fy(i,j,k) = zero
             else if (is_regular_cell(flag(i,j,k)) .or. apy(i,j,k).eq.one) then
                fy(i,j,k) = - dhy*bY(i,j,k)*(x(i,j,k) - x(i,j-1,k))
             else 
                fym = bY(i,j,k)*(x(i,j,k) - x(i,j-1,k))
                ii = i + int(sign(one,fcy(i,j,k,1)))
                kk = k + int(sign(one,fcy(i,j,k,2)))
                fracx = abs(fcy(i,j,k,1))*real(ior(ccm(ii,j-1,k),ccm(ii,j,k)),amrex_real)
                fracz = abs(fcy(i,j,k,2))*real(ior(ccm(i,j-1,kk),ccm(i,j,kk)),amrex_real)
                fym = (one-fracx)*(one-fracz)*fym + &
                     & fracx*(one-fracz)*bY(ii,j,k )*(x(ii,j,k )-x(ii,j-1,k )) + & 
                     & fracz*(one-fracx)*bY(i ,j,kk)*(x(i ,j,kk)-x(i ,j-1,kk)) + &
                     & fracx*     fracz *bY(ii,j,kk)*(x(ii,j,kk)-x(ii,j-1,kk))
                fy(i,j,k) = -dhy*fym
             endif
          enddo
       enddo
    enddo
    
    do       k = lo(3), hi(3)+1, kstride
       do    j = lo(2), hi(2) 
          do i = lo(1), hi(1) 
             if (apz(i,j,k) .eq. zero) then
                fz(i,j,k) = zero
             else if (is_regular_cell(flag(i,j,k)) .or. apz(i,j,k).eq.one) then
                fz(i,j,k) = - dhz*bZ(i,j,k)*(x(i,j,k) - x(i,j,k-1))
             else 
                fzm = bZ(i,j,k)*(x(i,j,k) - x(i,j,k-1))
                ii = i + int(sign(one,fcz(i,j,k,1)))
                jj = j + int(sign(one,fcz(i,j,k,2)))
                fracx = abs(fcz(i,j,k,1))*real(ior(ccm(ii,j,k-1),ccm(ii,j,k)),amrex_real)
                fracy = abs(fcz(i,j,k,2))*real(ior(ccm(i,jj,k-1),ccm(i,jj,k)),amrex_real)
                fzm = (one-fracx)*(one-fracy)*fzm + &
                     & fracx*(one-fracy)*bZ(ii,j ,k)*(x(ii,j ,k)-x(ii,j ,k-1)) + & 
                     & fracy*(one-fracx)*bZ(i ,jj,k)*(x(i ,jj,k)-x(i ,jj,k-1)) + &
                     & fracx*     fracy *bZ(ii,jj,k)*(x(ii,jj,k)-x(ii,jj,k-1))
                fz(i,j,k) = -dhz*fzm
             endif
          enddo
       enddo
    enddo
  end subroutine amrex_mlebabeclap_flux

  subroutine amrex_mlebabeclap_grad(xlo, xhi, ylo, yhi, zlo, zhi, sol, slo, shi,            &
                                    gx, gxlo, gxhi, gy, gylo, gyhi, gz, gzlo, gzhi,         &
                                    apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi,      &
                                    fcx, cxlo, cxhi, fcy, cylo, cyhi, fcz, czlo, czhi,      &
                                    ccm, cmlo, cmhi, flag, flo, fhi, dxinv) &
                                    bind(c, name='amrex_mlebabeclap_grad')
    
    integer, dimension(3), intent(in) :: xlo, xhi, slo, shi, ylo, yhi, zlo, zhi, &
         flo, fhi, axlo, axhi, aylo, ayhi, azlo, azhi,         &
         gxlo,gxhi,gylo, gyhi, gzlo, gzhi,  cxlo, cxhi, cylo,  &
         cyhi, czlo, czhi, cmlo, cmhi
    real(amrex_real), intent(in)     :: dxinv(3) 
    
    real(amrex_real), intent(inout) ::   gx(gxlo(1):gxhi(1),gxlo(2):gxhi(2),gxlo(3):gxhi(3))
    real(amrex_real), intent(inout) ::   gy(gylo(1):gyhi(1),gylo(2):gyhi(2),gylo(3):gyhi(3))
    real(amrex_real), intent(inout) ::   gz(gzlo(1):gzhi(1),gzlo(2):gzhi(2),gzlo(3):gzhi(3))
    real(amrex_real), intent(in   ) ::  sol( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    integer         , intent(in   ) ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3)) 
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3)) 
    real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)) 
    real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(amrex_real), intent(in   ) ::  apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
    real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2) 
    real(amrex_real), intent(in   ) ::  fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2) 
    
    integer  :: i, j, k, ii, jj, kk 
    real(amrex_real) :: dhx, dhy, dhz, fxm, fym, fzm, fracx, fracy, fracz
    
    dhx = dxinv(1) 
    dhy = dxinv(2)  
    dhz = dxinv(3)
    
    do       k = xlo(3), xhi(3) 
       do    j = xlo(2), xhi(2) 
          do i = xlo(1), xhi(1) 
             if(apx(i,i,j) .eq. zero) then
                gx(i,j,k) = zero
             else if (is_regular_cell(flag(i,j,k)) .or. apx(i,j,k).eq.one) then
                gx(i,j,k) = dhx*(sol(i,j,k) - sol(i-1,j,k)) 
             else 
                fxm = sol(i,j,k) - sol(i-1,j,k)
                jj = j + int(sign(one, fcx(i,j,k,1)))
                kk = k + int(sign(one, fcx(i,j,k,2)))
                fracy = abs(fcx(i,j,k,1))*real(ior(ccm(i-1,jj,k),ccm(i,jj,k)),amrex_real)
                fracz = abs(fcx(i,j,k,2))*real(ior(ccm(i-1,j,kk),ccm(i,j,kk)),amrex_real)
                fxm = (one-fracy)*(one-fracz)*fxm + &
                     & fracy*(one-fracz)*(sol(i,jj,k )-sol(i-1,jj,k )) + & 
                     & fracz*(one-fracy)*(sol(i,j ,kk)-sol(i-1,j ,kk)) + &
                     & fracy*     fracz *(sol(i,jj,kk)-sol(i-1,jj,kk))
                gx(i,j,k) = dhx*fxm
             endif
          enddo
       enddo
    enddo
    
    do       k = ylo(3), yhi(3) 
       do    j = ylo(2), yhi(2) 
          do i = ylo(1), yhi(1) 
             if(apy(i,j,k) .eq. zero) then
                gy(i,j,k) = zero
             else if (is_regular_cell(flag(i,j,k)) .or. apy(i,j,k).eq.one) then
                gy(i,j,k) = dhy*(sol(i,j,k) - sol(i,j-1,k))
             else 
                fym = (sol(i,j,k) - sol(i,j-1,k))
                ii = i + int(sign(one,fcy(i,j,k,1)))
                kk = k + int(sign(one,fcy(i,j,k,2)))
                fracx = abs(fcy(i,j,k,1))*real(ior(ccm(ii,j-1,k),ccm(ii,j,k)),amrex_real)
                fracz = abs(fcy(i,j,k,2))*real(ior(ccm(i,j-1,kk),ccm(i,j,kk)),amrex_real)
                fym = (one-fracx)*(one-fracz)*fym + &
                     & fracx*(one-fracz)*(sol(ii,j,k )-sol(ii,j-1,k )) + & 
                     & fracz*(one-fracx)*(sol(i ,j,kk)-sol(i ,j-1,kk)) + &
                     & fracx*     fracz *(sol(ii,j,kk)-sol(ii,j-1,kk))
                gy(i,j,k) = dhy*fym
             endif
          enddo
       enddo
    enddo
    
    do       k = zlo(3), zhi(3) 
       do    j = zlo(2), zhi(2) 
          do i = zlo(1), zhi(1) 
             if(apz(i,j,k) .eq. zero) then
                gz(i,j,k) = zero
             else if (is_regular_cell(flag(i,j,k)) .or. apz(i,j,k).eq.one) then
                gz(i,j,k) = dhz*(sol(i,j,k) - sol(i,j,k-1))
             else 
                fzm = (sol(i,j,k) - sol(i,j,k-1))
                ii = i + int(sign(one,fcz(i,j,k,1)))
                jj = j + int(sign(one,fcz(i,j,k,2)))
                fracx = abs(fcz(i,j,k,1))*real(ior(ccm(ii,j,k-1),ccm(ii,j,k)),amrex_real)
                fracy = abs(fcz(i,j,k,2))*real(ior(ccm(i,jj,k-1),ccm(i,jj,k)),amrex_real)
                fzm = (one-fracx)*(one-fracy)*fzm + &
                     & fracx*(one-fracy)*(sol(ii,j ,k)-sol(ii,j ,k-1)) + & 
                     & fracy*(one-fracx)*(sol(i ,jj,k)-sol(i ,jj,k-1)) + &
                     & fracx*     fracy *(sol(ii,jj,k)-sol(ii,jj,k-1))
                gz(i,j,k) = dhz*fzm
             endif
          enddo
       enddo
    enddo
  end subroutine amrex_mlebabeclap_grad

  subroutine compute_dphidn_3d (dphidn, dxinv, i, j, k, &
        phi,  p_lo, p_hi,     &
        flag,  flo,  fhi,     &
        bct, phib, anrmx, anrmy, anrmz, vf)

      ! Cell indices 
      integer, intent(in   ) :: i, j, k

      ! Grid spacing
      real(amrex_real),       intent(in   ) :: dxinv(3)

      ! Array bounds
      integer, intent(in   ) :: p_lo(3), p_hi(3)
      integer, intent(in   ) ::  flo(3),  fhi(3)

      ! Arrays
      real(amrex_real),  intent(in   ) ::                            &
           & phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
      integer, intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3)) 

      real(amrex_real),  intent(in   ) :: bct(3), phib
      real(amrex_real),  intent(in   ) :: anrmx, anrmy, anrmz, vf

      real(amrex_real),        intent(  out) :: dphidn

      ! Local variable
      real(amrex_real) :: bctx, bcty, bctz
      real(amrex_real) :: phig, gx, gy, gz, dg, dx_eb, gxy, gxz, gyz, gxyz
      real(amrex_real) :: sx, sy, sz
      integer          :: ii, jj, kk

      bctx = bct(1)
      bcty = bct(2)
      bctz = bct(3)

      dx_eb = amrex_get_dx_eb(vf)

      dg = dx_eb / max(abs(anrmx),abs(anrmy),abs(anrmz))
      gx = bctx - dg*anrmx
      gy = bcty - dg*anrmy
      gz = bctz - dg*anrmz
      sx =  sign(one,anrmx)
      sy =  sign(one,anrmy)
      sz =  sign(one,anrmz)
      ii = i - int(sx)
      jj = j - int(sy)
      kk = k - int(sz)

      gx = sx*gx
      gy = sy*gy
      gz = sz*gz
      gxy = gx*gy
      gxz = gx*gz
      gyz = gy*gz
      gxyz = gx*gy*gz
      phig = (one+gx+gy+gz+gxy+gxz+gyz+gxyz) * phi(i ,j ,k ) &
           + (-gz - gxz - gyz - gxyz)        * phi(i ,j ,kk) &
           + (-gy - gxy - gyz - gxyz)        * phi(i ,jj,k ) &
           + (gyz + gxyz)                    * phi(i ,jj,kk) &
           + (-gx - gxy - gxz - gxyz)        * phi(ii,j ,k ) &
           + (gxz + gxyz)                    * phi(ii,j ,kk) &
           + (gxy + gxyz)                    * phi(ii,jj,k ) &
           + (-gxyz)                         * phi(ii,jj,kk)

      dphidn = (phib-phig)/dg

  end subroutine compute_dphidn_3d

  subroutine compute_dphidn_3d_ho (dphidn, dxinv, i, j, k, &
        phi,  p_lo, p_hi,     &
        flag,  flo,  fhi,     &
        bct, phib, anrmx, anrmy, anrmz)

      ! Cell indices 
      integer, intent(in   ) :: i, j, k

      ! Grid spacing
      real(amrex_real),       intent(in   ) :: dxinv(3)

      ! Array bounds
      integer, intent(in   ) :: p_lo(3), p_hi(3)
      integer, intent(in   ) ::  flo(3),  fhi(3)

      ! Arrays
      real(amrex_real),  intent(in   ) ::                            &
           & phi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))

      integer, intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3)) 

      real(amrex_real),  intent(in   ) :: bct(3), phib
      real(amrex_real),  intent(in   ) :: anrmx, anrmy, anrmz

      real(amrex_real),        intent(  out) :: dphidn

      ! Local variable
      real(amrex_real) :: anrm
      real(amrex_real) :: xit, yit, zit, s, s2
      real(amrex_real) :: d1, d2, ddinv
      real(amrex_real) :: u1, u2
      integer          :: ixit, iyit, izit, is, is2, ivar

      if (abs(anrmx).ge.abs(anrmy) .and. abs(anrmx).ge.abs(anrmz)) then
         anrm = anrmx
         ivar = 1
      else if (abs(anrmy).ge.abs(anrmx) .and. abs(anrmy).ge.abs(anrmz)) then
         anrm = anrmy
         ivar = 2
      else 
         anrm = anrmz
         ivar = 3
      end if

      ! s is -1. or 1.
        s = sign(one,-anrm)
       s2 = 2*s

      ! is is -1 or 1
      is  = nint(s)
      is2 = 2*is

      d1 = (bct(ivar) - s ) * (one/anrm)
      d2 = (bct(ivar) - s2) * (one/anrm)

      if (abs(anrmx).ge.abs(anrmy) .and. abs(anrmx).ge.abs(anrmz)) then
         ! y-z plane: x = const
         ! the equation for the line:  x = bct(1) - d*anrmx
         !                             y = bct(2) - d*anrmy
         !                             z = bct(3) - d*anrmz

         !
         ! the line intersects the y-z plane (x = s) at ...
         !
         yit = bct(2) - d1*anrmy
         zit = bct(3) - d1*anrmz
         iyit = j + nint(yit)
         izit = k + nint(zit)
         yit = yit - nint(yit)  ! shift so that the center of the nine cells are (0.,0.)
         zit = zit - nint(zit)

         call interp2d(u1,yit,zit,phi(i+is,iyit-1:iyit+1,izit-1:izit+1),flag(i+is,iyit-1:iyit+1,izit-1:izit+1))

         !
         ! the line intersects the y-z plane (x = 2*s) at ...
         !
         yit = bct(2) - d2*anrmy
         zit = bct(3) - d2*anrmz
         iyit = j + nint(yit)
         izit = k + nint(zit)
         yit = yit - nint(yit)  ! shift so that the center of the nine cells are (0.,0.)
         zit = zit - nint(zit)

         call interp2d(u2,yit,zit,phi(i+is2,iyit-1:iyit+1,izit-1:izit+1),flag(i+is2,iyit-1:iyit+1,izit-1:izit+1))

      else if (abs(anrmy).ge.abs(anrmx) .and. abs(anrmy).ge.abs(anrmz)) then

         xit = bct(1) - d1*anrmx
         zit = bct(3) - d1*anrmz
         ixit = i + nint(xit)
         izit = k + nint(zit)
         xit = xit - nint(xit)
         zit = zit - nint(zit)

         call interp2d(u1,xit,zit,phi(ixit-1:ixit+1,j+is,izit-1:izit+1),flag(ixit-1:ixit+1,j+is,izit-1:izit+1))

         xit = bct(1) - d2*anrmx
         zit = bct(3) - d2*anrmz
         ixit = i + nint(xit)
         izit = k + nint(zit)
         xit = xit - nint(xit)
         zit = zit - nint(zit)

         call interp2d(u2,xit,zit,phi(ixit-1:ixit+1,j+is2,izit-1:izit+1),flag(ixit-1:ixit+1,j+is2,izit-1:izit+1))

      else

         xit = bct(1) - d1*anrmx
         yit = bct(2) - d1*anrmy
         ixit = i + nint(xit)
         iyit = j + nint(yit)
         xit = xit - nint(xit)
         yit = yit - nint(yit)

         call interp2d(u1,xit,yit,phi(ixit-1:ixit+1,iyit-1:iyit+1,k+is),flag(ixit-1:ixit+1,iyit-1:iyit+1,k+is))

         xit = bct(1) - d2*anrmx
         yit = bct(2) - d2*anrmy
         ixit = i + nint(xit)
         iyit = j + nint(yit)
         xit = xit - nint(xit)
         yit = yit - nint(yit)

         call interp2d(u2,xit,yit,phi(ixit-1:ixit+1,iyit-1:iyit+1,k+is2),flag(ixit-1:ixit+1,iyit-1:iyit+1,k+is2))

      end if

      !
      ! compute derivatives on the wall given that phi is zero on the wall.
      !
      ddinv = one/(d1*d2*(d2-d1))
      dphidn = -ddinv*( d2*d2*(u1-phib) - d1*d1*(u2-phib) )  ! note that the normal vector points toward the wall

  end subroutine compute_dphidn_3d_ho

  subroutine interp2d(phi_interp,yit,zit,v,flag)

      real(amrex_real), intent(in   ) :: yit,zit,v(3,3)
      integer         , intent(in   ) :: flag(3,3)
      real(amrex_real), intent(  out) :: phi_interp

      real(amrex_real)             :: cym,cy0,cyp,czm,cz0,czp
      real(amrex_real)             :: val(3)

      ! Coefficents for quadratic interpolation
      cym = half*yit*(yit-one)
      cy0 = one-yit*yit
      cyp = half*yit*(yit+one)

      ! Coefficents for quadratic interpolation
      czm = half*zit*(zit-one)
      cz0 = one-zit*zit
      czp = half*zit*(zit+one)

      if (any(is_covered_cell(flag(1,1:3)))) then

         val(1:3) = (one-yit)*v(2,1:3) + yit*v(3,1:3)

         if (any(is_covered_cell(flag(2:3,1)))) then
            phi_interp = (one-yit)*val(2) + yit*val(3)
         else if (any(is_covered_cell(flag(2:3,3)))) then
            phi_interp = -yit*val(1) + (one+yit)*val(2)
         else 
            phi_interp =  czm*val(1) + cz0*val(2) + czp*val(3)
         end if

      else if (any(is_covered_cell(flag(3,1:3)))) then

         val(1:3) = -yit*v(1,1:3) + (one+yit)*v(2,1:3)

         if (any(is_covered_cell(flag(1:2,1)))) then
            phi_interp = (one-yit)*val(2) + yit*val(3)
         else if (any(is_covered_cell(flag(1:2,3)))) then
            phi_interp = -yit*val(1) + (one+yit)*val(2)
         else 
            phi_interp =  czm*val(1) + cz0*val(2) + czp*val(3)
         end if

      else

         val(1:3) =  cym*v(1,1:3) + cy0*v(2,1:3) + cyp*v(3,1:3)

         if (any(is_covered_cell(flag(1:3,1)))) then
            phi_interp = (one-yit)*val(2) + yit*val(3)
         else if (any(is_covered_cell(flag(1:3,3)))) then
            phi_interp = -yit*val(1) + (one+yit)*val(2)
         else 
            phi_interp =  czm*val(1) + cz0*val(2) + czp*val(3)
         end if

      end if

  end subroutine interp2d

end module amrex_mlebabeclap_3d_module
