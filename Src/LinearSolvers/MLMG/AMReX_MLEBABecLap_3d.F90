
module amrex_mlebabeclap_3d_module

  use amrex_error_module 
  use amrex_constants_module, only: zero, one 
  use amrex_fort_module, only : amrex_real
  use amrex_ebcellflag_module, only: is_regular_cell, is_covered_cell, is_single_valued_cell, &
       get_neighbor_cells_int_single
  implicit none

  private
  public :: amrex_mlebabeclap_adotx, amrex_mlebabeclap_gsrb, amrex_mlebabeclap_normalize, & 
       amrex_eb_mg_interp

contains

  subroutine amrex_mlebabeclap_adotx(lo, hi, y, ylo, yhi, x, xlo, xhi, & 
       a, alo, ahi, bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, flag, flo, fhi, & 
       vfrc, vlo, vhi, apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi, fcx, cxlo, cxhi, &
       fcy, cylo, cyhi, fcz, czlo, czhi, cntr, ctlo, cthi, dxinv, alpha, beta) & 
       bind(c, name='amrex_mlebabeclap_adotx') 
   integer, dimension(3), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, alo, ahi, bxlo,&
      bxhi, bylo, byhi, bzlo, bzhi, flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, azlo, azhi, & 
      cxlo, cxhi, cylo, cyhi, czlo, czhi, ctlo, cthi
   real(amrex_real), intent(in) :: dxinv(3) 
   real(amrex_real), value, intent(in) :: alpha, beta
   real(amrex_real), intent(inout) ::    y( ylo(1): yhi(1), ylo(2): yhi(2), ylo(3): yhi(3))
   real(amrex_real), intent(in   ) ::    x( xlo(1): xhi(1), xlo(2): xhi(2), xlo(3): xhi(3))
   real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
   real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
   real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
   real(amrex_real), intent(in   ) ::   bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
   integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3)) 
   real(amrex_real), intent(in   ) :: vfrc( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3)) 
   real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)) 
   real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
   real(amrex_real), intent(in   ) ::  apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
   real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
   real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2) 
   real(amrex_real), intent(in   ) ::  fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2) 
   real(amrex_real), intent(in   ) :: cntr(ctlo(1):cthi(1),ctlo(2):cthi(2),ctlo(3):cthi(3),3)
   integer  :: i, j, k, ii, jj, kk 
   real(amrex_real) :: dhx, dhy, dhz, fxm, fxp, fym, fyp, fzm, fzp, fracx, fracy, fracz, cx, cy, cz 
   real(amrex_real) :: face_cen(6,2), dxa, dya, dza   

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
                    fracy = abs(fcx(i,j,k,1))
                    fracz = abs(fcx(i,j,k,2))
                    jj = j + int(sign(one, fcx(i,j,k,1)))
                    kk = k + int(sign(one, fcx(i,j,k,2)))
                    fxm = (one-fracy-fracz)*fxm + fracy*bX(i,jj,k)*(x(i,jj,k)-x(i-1,jj,k)) + & 
                           fracz*bX(i,j,kk)*(x(i,j,kk)-x(i-1,j,kk))
                endif 

                fxp = bX(i+1,j,k)*(x(i+1,j,k) - x(i,j,k))
                if (apx(i+1,j,k).ne.zero.and.apx(i+1,j,k).ne.one) then 
                    fracy = abs(fcx(i+1,j,k,1))
                    fracz = abs(fcx(i+1,j,k,2))
                    jj = j + int(sign(one,fcx(i+1,j,k,1)))
                    kk = k + int(sign(one,fcx(i+1,j,k,2)))
                    fxp = (one-fracy-fracz)*fxp + fracy*bX(i+1,jj,k)*(x(i+1,jj,k)-x(i,jj,k)) + & 
                           fracz*bX(i+1,j,kk)*(x(i+1,j,kk)-x(i,j,kk))
                endif 

                fym = bY(i,j,k)*(x(i,j,k) - x(i,j-1,k))
                if (apy(i,j,k).ne.zero.and.apy(i,j,k).ne.one) then 
                    fracx = abs(fcy(i,j,k,1))
                    fracz = abs(fcy(i,j,k,2))
                    ii = i + int(sign(one,fcy(i,j,k,1)))
                    kk = k + int(sign(one,fcy(i,j,k,2)))
                    fym = (one-fracx-fracz)*fym + fracx*bY(ii,j,k)*(x(ii,j,k) - x(ii,j-1,k)) + & 
                           fracz*bY(ii,j,kk)*(x(i,j,kk)-x(i,j-1,kk))
                endif 

                fyp = bY(i,j+1,k)*(x(i,j+1,k) - x(i,j,k))
                if (apy(i,j+1,k).ne.zero.and.apy(i,j+1,k).ne.one) then 
                    fracx = abs(fcy(i,j+1,k,1))
                    fracz = abs(fcy(i,j+1,k,2))
                    ii = i + int(sign(one,fcy(i,j+1,k,1)))
                    kk = k + int(sign(one,fcy(i,j+1,k,2)))
                    fyp = (one-fracx-fracz)*fyp + fracx*bY(ii,j+1,k)*(x(ii,j+1,k) - x(ii,j,k)) + &
                           fracz*bY(i,j+1,kk)*(x(i,j+1,kk)-x(i,j,kk))
                endif 

                fzm = bZ(i,j,k)*(x(i,j,k) - x(i,j,k-1))
                if (apz(i,j,k).ne.zero.and.apz(i,j,k).ne.one) then 
                    fracx = abs(fcz(i,j,k,1))
                    fracy = abs(fcz(i,j,k,2))
                    ii = i + int(sign(one,fcz(i,j,k,1)))
                    jj = j + int(sign(one,fcz(i,j,k,2)))
                    fzm = (one-fracx-fracy)*fzm + fracx*bZ(ii,jj,k)*(x(ii,j,k) - x(ii,j,k-1)) + & 
                           fracy*bZ(i,jj,k)*(x(i,jj,k)-x(i,jj,k-1))
                endif 

                fzp = bZ(i,j,k+1)*(x(i,j,k+1) - x(i,j,k))
                if (apz(i,j,k+1).ne.zero.and.apz(i,j,k+1).ne.one) then 
                    fracx = abs(fcz(i,j,k+1,1))
                    fracy = abs(fcz(i,j,k+1,2))
                    ii = i + int(sign(one,fcz(i,j,k+1,1)))
                    jj = j + int(sign(one,fcz(i,j,k+1,2)))
                    fzp = (one-fracx-fracy)*fzp + fracx*bZ(ii,j,k+1)*(x(ii,j,k+1) - x(ii,j,k)) + &
                           fracy*bZ(i,jj,k+1)*(x(i,jj,k+1)-x(i,jj,k))
                endif 

                y(i,j,k) = (one/vfrc(i,j,k))*&
                       (dhx*(apx(i,j,k)*fxm - apx(i+1,j,k)*fxp) + &
                        dhy*(apy(i,j,k)*fym - apy(i,j+1,k)*fyp) + &
                        dhz*(apz(i,j,k)*fzm - apz(i,j,k+1)*fzp))

               if(alpha .ne. zero) then
                dxa = zero 
                dya = zero 
                dza = zero 
                cx  = cntr(i,j,k,1) 
                cy  = cntr(i,j,k,2) 
                cz  = cntr(i,j,k,3) 

                !derivatives of "ax"
                if(cx.ge.zero) then 
                   dxa = dxinv(1)*(a(i+1,j,k)*x(i+1,j,k) - a(i,j,k)*x(i,j,k))
                else
                   dxa = dxinv(1)*(a(i,j,k)*x(i,j,k) - a(i-1,j,k)*x(i-1,j,k)) 
                endif
                if(cy.ge.zero) then 
                   dya = dxinv(2)*(a(i,j+1,k)*x(i,j+1,k) - a(i,j,k)*x(i,j,k)) 
                else
                   dya = dxinv(2)*(a(i,j,k)*x(i,j,k) - a(i,j-1,k)*x(i,j-1,k)) 
                endif
                if(cz.ge.zero) then 
                   dza = dxinv(3)*(a(i,j,k+1)*x(i,j,k+1) - a(i,j,k)*x(i,j,k)) 
                else 
                   dza = dxinv(3)*(a(i,j,k)*x(i,j,k) - a(i,j,k-1)*x(i,j,k-1)) 
                endif 
       
                y(i,j,k) = y(i,j,k) + alpha*(a(i,j,k) + cx*dxa + cy*dya + cz*dza)*x(i,j,k)
               endif
              endif 
          enddo
       enddo 
   enddo 
  end subroutine amrex_mlebabeclap_adotx

  subroutine amrex_mlebabeclap_gsrb(lo, hi, phi, hlo, hhi, rhs, rlo, rhi, a, alo, ahi, & 
     bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, & 
     m0, m0lo, m0hi, m2, m2lo, m2hi, m4, m4lo, m4hi, & 
     m1, m1lo, m1hi, m3, m3lo, m3hi, m5, m5lo, m5hi, & 
     f0, f0lo, f0hi, f2, f2lo, f2hi, f4, f4lo, f4hi, & 
     f1, f1lo, f1hi, f3, f3lo, f3hi, f5, f5lo, f5hi, & 
     flag, flo, fhi, vfrc, vlo, vhi, & 
     apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
     fcz, czlo, czhi, dxinv, alpha, beta, redblack) & 
     bind(c,name='amrex_mlebabeclap_gsrb') 

    integer, dimension(3), intent(in) :: lo, hi, hlo, hhi, rlo, rhi, alo, ahi, bxlo, bxhi, bylo, byhi, &
         bzlo, bzhi, m0lo, m0hi, m1lo, m1hi, m2lo, m2hi, m3lo, m3hi, m4lo, m4hi, m5lo, m5hi,  &
         f0lo, f0hi, f1lo, f1hi, f2lo, f2hi, f3lo, f3hi, f4lo, f4hi ,f5lo, f5hi, flo, fhi, vlo, vhi,&
         axlo, axhi, aylo, ayhi, azlo, azhi, cxlo, cxhi, cylo, cyhi, czlo, czhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), value, intent(in) :: alpha, beta
    integer         , value, intent(in) :: redblack
    real(amrex_real), intent(inout) ::  phi( hlo(1): hhi(1), hlo(2): hhi(2), hlo(3): hhi(3)  )
    real(amrex_real), intent(in   ) ::  rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3)  )
    real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3)  )
    real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3)  )
    real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3)  )
    real(amrex_real), intent(in   ) ::   bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3)  )
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

    integer :: i, j, k, ioff, ii, jj, kk
    real(amrex_real) :: cf0, cf1, cf2, cf3, cf4, cf5,  delta, gamma, rho 
    real(amrex_real) :: dhx, dhy, dhz, fxm, fxp, fym, fyp, fzm, fzp, fracx, fracy, fracz
    real(amrex_real) :: sxm, sxp, sym, syp, szm, szp
    real(amrex_real) :: face_cen(6,2), dxa, dya, dza, cx, cy, cz 
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
                delta = dhx*(bX(i,j,k)*cf0 + bX(i+1,j,k)*cf3) & 
                      + dhy*(bY(i,j,k)*cf1 + bY(i,j+1,k)*cf4) & 
                      + dhz*(bZ(i,j,k)*cf2 + bZ(i,j,k+1)*cf5)

                if(is_regular_cell(flag(i,j,k))) then 
                  
                   gamma = alpha*a(i,j,k) & 
                         + dhx*(bX(i+1,j,k) + bX(i,j,k)) & 
                         + dhy*(bY(i,j+1,k) + bY(i,j,k)) &
                         + dhz*(bZ(i,j,k+1) + bZ(i,j,k)) 

                   rho   = dhx*(bX(i+1,j,k)*phi(i+1,j,k) + bX(i,j,k)*phi(i-1,j,k)) & 
                         + dhy*(bY(i,j+1,k)*phi(i,j+1,k) + bY(i,j,k)*phi(i,j-1,k)) & 
                         + dhz*(bZ(i,j,k+1)*phi(i,j,k+1) + bZ(i,j,k)*phi(i,j,k-1)) 
                else 
                   dxa = zero
                   dya = zero 
                   dza = zero 
 
                   fxm = -bX(i,j,k)*phi(i-1,j,k) 
                   sxm =  bX(i,j,k) 
                   if(apx(i,j,k).ne.zero .and. apx(i,j,k).ne.one) then 
                       fracy = abs(fcx(i,j,k,1))
                       fracz = abs(fcx(i,j,k,2)) 
                       jj = j + int(sign(one, fcx(i,j,k,1)))
                       kk = k + int(sign(one, fcx(i,j,k,2))) 
                       fxm = (one-fracy-fracz)*fxm + fracy*bX(i,jj,k)*(phi(i,jj,k)-phi(i-1,jj,k)) &
                           +  fracz*bX(i,j,kk)*(phi(i,j,kk)-phi(i-1,j,kk))
                       sxm = (one-fracy-fracz)*sxm
                   end if
                   
                   fxp =  bX(i+1,j,k)*phi(i+1,j,k) 
                   sxp = -bX(i+1,j,k)
                   if(apx(i+1,j,k).ne.zero.and.apx(i+1,j,k).ne.one) then 
                       fracy = abs(fcx(i+1,j,k,1)) 
                       fracz = abs(fcx(i+1,j,k,2)) 
                       jj = j + int(sign(one, fcx(i+1,j,k,1)))
                       kk = k + int(sign(one, fcx(i+1,j,k,2)))
                       fxp = (one-fracy-fracz)*fxp + fracy*bX(i+1,jj,k)*(phi(i+1,jj,k)-phi(i,jj,k)) & 
                           + fracz*bX(i+1,j,kk)*(phi(i+1,j,kk)-phi(i,j,kk))
                       sxp = (one-fracy-fracz)*sxp
                   end if 
                   
                   fym = -bY(i,j,k)*phi(i,j-1,k)
                   sym =  bY(i,j,k)
                   if(apy(i,j,k).ne.zero.and.apy(i,j,k).ne.one) then 
                      fracx = abs(fcy(i,j,k,1))
                      fracz = abs(fcy(i,j,k,2))
                      ii = i + int(sign(one,fcy(i,j,k,1)))
                      kk = k + int(sign(one,fcy(i,j,k,2))) 
                      fym = (one-fracx-fracz)*fym + fracx*bY(ii,j,k)*(phi(ii,j,k)-phi(ii,j-1,k)) & 
                          +  fracz*bY(i,j,kk)*(phi(i,j,kk)-phi(i,j-1,kk))
                      sym = (one-fracx-fracz)*sym
                   endif
 
                   fyp =  bY(i,j+1,k)*phi(i,j+1,k)
                   syp = -bY(i,j+1,k)
                   if(apy(i,j+1,k).ne.zero.and.apy(i,j+1,k).ne.one) then 
                      fracx = abs(fcy(i,j+1,k,1))
                      fracz = abs(fcy(i,j+1,k,2)) 
                      ii = i + int(sign(one,fcy(i,j+1,k,1)))
                      kk = k + int(sign(one,fcy(i,j+1,k,2)))
                      fyp = (one-fracx-fracz)*fyp + fracx*bY(ii,j+1,k)*(phi(ii,j+1,k)-phi(ii,j,k)) &
                          +  fracz*bY(i,j+1,kk)*(phi(i,j+1,kk)-phi(i,j,kk))
                      syp = (one-fracx-fracz)*syp
                   end if 
 
                   fzm = -bZ(i,j,k)*phi(i,j,k-1)
                   szm =  bZ(i,j,k)
                   if(apz(i,j,k).ne.zero.and.apz(i,j,k).ne.one) then 
                      fracx = abs(fcz(i,j,k,1))
                      fracy = abs(fcz(i,j,k,2)) 
                      ii = i + int(sign(one,fcz(i,j,k,1)))
                      jj = j + int(sign(one,fcz(i,j,k,2)))
                      fzm = (one-fracx-fracy)*fzm + fracx*bZ(ii,j,k)*(phi(ii,j,k)-phi(ii,j,k-1)) & 
                          +  fracy*bZ(i,jj,k)*(phi(i,jj,k)-phi(i,jj,k-1))
                      szm = (one-fracx-fracy)*szm
                    endif
               
                    fzp =  bZ(i,j,k+1)*phi(i,j,k+1) 
                    szp = -bZ(i,j,k+1)
                    if(apz(i,j,k+1).ne.zero.and.apz(i,j,k+1).ne.one) then 
                       fracx = abs(fcz(i,j,k+1,1))
                       fracy = abs(fcz(i,j,k+1,2))
                       ii = i + int(sign(one,fcz(i,j,k+1,1)))
                       jj = j + int(sign(one,fcz(i,j,k+1,2)))
                       fzp = (one-fracx-fracy)*fzp + fracx*bZ(ii,j,k+1)*(phi(ii,j,k+1)-phi(ii,j,k)) &
                           +  fracy*bZ(i,jj,k+1)*(phi(i,jj,k+1)-phi(i,jj,k))
                       szp = (one-fracx-fracy)*szp
                    end if 

                  
                    gamma = alpha*a(i,j,k) + & 
                           (one/vfrc(i,j,k))*(dhx*(apx(i,j,k)*sxm-apx(i+1,j,k)*sxp) + &
                            dhy*(apy(i,j,k)*sym-apy(i,j+1,k)*syp) + &
                            dhz*(apz(i,j,k)*szm-apz(i,j,k+1)*szp))

                    rho = -(one/vfrc(i,j,k)) * & 
                           (dhx*(apx(i,j,k)*fxm-apx(i+1,j,k)*fxp) + &
                            dhy*(apy(i,j,k)*fym-apy(i,j+1,k)*fyp) + &
                            dhz*(apz(i,j,k)*fzm-apz(i,j,k+1)*fzp))
                  end if 

                  phi(i,j,k) = (rhs(i,j,k) + rho - phi(i,j,k)*delta)/(gamma -delta)
              endif
          end do 
       end do 
    enddo        
 end subroutine amrex_mlebabeclap_gsrb


  subroutine amrex_mlebabeclap_normalize (lo, hi, x, xlo, xhi, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, flag, flo, fhi, vfrc, vlo, vhi, &
       apx, axlo, axhi, apy, aylo, ayhi,apz, azlo, azhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
       fcz, czlo, czhi, dxinv, alpha, beta) &
       bind(c,name='amrex_mlebabeclap_normalize')

    integer, dimension(3), intent(in) :: lo, hi, xlo, xhi, alo, ahi, bxlo, bxhi, bylo, byhi, bzlo, bzhi, &
              flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi,azlo, azhi, cxlo, cxhi, cylo, cyhi, czlo, czhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), value, intent(in) :: alpha, beta
    real(amrex_real), intent(inout) ::    x( xlo(1): xhi(1), xlo(2): xhi(2), xlo(3): xhi(3)  )
    real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3)  )
    real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bzhi(3)  )
    real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3)  )
    real(amrex_real), intent(in   ) ::   bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3)  )
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3)  )
    real(amrex_real), intent(in   ) :: vfrc( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3)  )
    real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)  )
    real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)  )
    real(amrex_real), intent(in   ) ::  apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3)  )
    real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
    real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2)
    real(amrex_real), intent(in   ) ::  fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2)

    integer :: i, j, k
    real(amrex_real) :: dhx, dhy, dhz, sxm, sxp, sym, syp, szm, szp, gamma
    real(amrex_real) :: face_cen(6,2), dxa, dya, dza, cx, cy, cz 

    dhx = beta*dxinv(1)*dxinv(1)
    dhy = beta*dxinv(2)*dxinv(2)
    dhz = beta*dxinv(3)*dxinv(3) 

    do      k = lo(3), hi(3)
      do    j = lo(2), hi(2)
        do  i = lo(1), hi(1)
          if (is_regular_cell(flag(i,j,k))) then
             x(i,j,k) = x(i,j,k) / (alpha*a(i,j,k) + dhx*(bX(i,j,k)+bX(i+1,j,k))   &
                      + dhy*(bY(i,j,k)+bY(i,j+1,k)) + dhz*(bZ(i,j,k)+bZ(i,j,k+1)))
          else if (is_single_valued_cell(flag(i,j,k))) then
             sxm =  bX(i,j,k  ) * (one-(abs(fcx(i  ,j,k,1))+abs(fcx(i  ,j,k,2))))
             sxp = -bX(i+1,j,k) * (one-(abs(fcx(i+1,j,k,1))+abs(fcx(i+1,j,k,2))))
             sym =  bY(i,j,k  ) * (one-(abs(fcy(i,j  ,k,1))+abs(fcy(i,j  ,k,2))))
             syp = -bY(i,j+1,k) * (one-(abs(fcy(i,j+1,k,1))+abs(fcy(i,j+1,k,2))))
             szm =  bZ(i,j,k  ) * (one-(abs(fcz(i,j,k  ,1))+abs(fcz(i,j,k  ,2))))
             szp = -bZ(i,j,k+1) * (one-(abs(fcz(i,j,k+1,1))+abs(fcz(i,j,k+1,2))))

             gamma =  alpha*a(i,j,k) + (one/vfrc(i,j,k)) * &
                  (dhx*(apx(i,j,k)*sxm-apx(i+1,j,k)*sxp) + &
                   dhy*(apy(i,j,k)*sym-apy(i,j+1,k)*syp) + & 
                   dhz*(apz(i,j,k)*szm-apz(i,j,k+1)*szp))

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


end module amrex_mlebabeclap_3d_module
