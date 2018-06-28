
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
       vfrc, vlo, vhi, apx, axlo, axhi, apy, aylo, ayhi, az, azlo, azhi, fcx, cxlo, cxhi, &
       fcy, cylo, cyhi, fcz, czlo, czhi, dxinv, alpha, beta) & 
       bind(c, name='amrex_mlebabeclap_adotx') 
   integer, dimension(3), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, alo, ahi, bxlo,&
      bxhi, bylo, byhi, bzlo, bzhi, flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, azlo, azhi, & 
      cxlo, cxhi, cylo, cyhi, czlo, czhi
   real(amrex_real), intent(in) :: dxinv(3) 
   real(amrex_real), value, intent(in) :: alpha, beta
   real(amrex_real), intent(inout) ::    y( ylo(1): yhi(1), ylo(2): yhi(2), ylo(3): yhi(3))
   real(amrex_real), intent(in   ) ::    x( xlo(1): xhi(1), xlo(2): xhi(2), xlo(3): xhi(3))
   real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
   real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
   real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
   real(amrex_real), intent(in   ) ::   bz(bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
   integer           intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3)) 
   real(amrex_real), intent(in   ) :: vfrc( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3)) 
   real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)) 
   real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
   real(amrex_real), intent(in   ) ::  apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
   real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3))
   real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3)) 
   real(amrex_real), intent(in   ) ::  fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3)) 

   integer  :: i, j, k, ii, jj, kk 
   real(amrex_real) :: dhx, dhy, dhz, fxm, fxp, fym, fyp, fzm, fzp, fracx, fracy, fracz 
   logical  :: xm, xp, ym, yp, zm, zp

   dhx = beta*dxinv(1)*dxinv(1) 
   dhy = beta*dxinv(2)*dxinv(2)  
   dhz = beta*dxinv(3)*dxinv(3) 

   do         k = lo(3), hi(3) 
       do     j = lo(2), hi(2) 
           do i = lo(1), hi(1) 
              xm = apx(i,j,k).ne.zero.and.apx(i,j,k).ne.one
              xp = apx(i+1,j,k).ne.zero.and.apx(i+1,j,k).ne.one 
              ym = apy(i,j,k).ne.zero.and.apy(i,j,k).ne.one
              yp = apy(i,j+1,k).ne.zero.and.apy(i,j+1,k).ne.one 
              zm = apz(i,j,k).ne.zero.and.apz(i,j,k).ne.one
              zp = apz(i,j,k+1).ne.zero.and.apz(i,j,k+1).ne.one 
             if(is_covered_cell(flag(i,j,k))) then 
                 y(i,j,k) = zero
              else if (is_regular_cell(flag(i,j,k))) then 
                 y(i,j,k) = alpha*a(i,j,k)*x(i,j,k) & 
                            - dhx * (bX(i+1,j,k)*(x(i+1,j,k) - x(i  ,j,k))  & 
                                    -bX(i  ,j,k)*(x(i  ,j,k) - x(i-1,j,k))) & 
                            - dhy * (bY(i,j+1,k)*(x(i,j+1,k) - x(i,j  ,k))  &
                                    -bY(i,j  ,k)*(x(i,j  ,k) - x(i,j-1,k))) &
                            - dhz * (bZ(i,j,k+1)*(x(i,j,k+1) - x(i,j,k  ))  & 
                                    -bZ(i,j,k  )*(x(i,j,k  ) - x(i,j,k-1)))
              else 
                fxm = bX(i,j,k)*(x(i,j,k) - x(i-1,j,k))
                if (xm) then 
                    fracy = abs(fcx(i,j,k))
                    jj = j
                    kk = k 
                    if(ym) kk = k + int(sign(one, fcx(i,j,k)))
                    if(zm) jj = j + int(sign(one, fcx(i,j,k)))
                   fxm = (one-fracy)*fxm + fracy*bX(i,jj,kk)*(x(i,jj,kk)-x(i-1,jj,kk))
                endif 

                fxp = bX(i+1,j,k)*(x(i+1,j,k) - x(i,j,k))
                if (xp) then 
                    fracy = abs(fcx(i+1,j,k))
                    jj = j 
                    kk = k 
                    if(yp) kk = k + int(sign(one,fcx(i+1,j,k)))
                    if(zp) jj = j + int(sign(one,fcx(i+1,j,k)))
                    fxp = (one-fracy)*fxp + fracy*bX(i+1,jj,kk)*(x(i+1,jj,kk)-x(i,jj,kk))
                endif 

                fym = bY(i,j,k)*(x(i,j,k) - x(i,j-1,k))
                if (ym) then 
                    fracz = abs(fcy(i,j,k))
                    ii = i 
                    kk = k 
                    if(zm) ii = i + int(sign(one,fcy(i,j,k)))
                    if(xm) kk = k + int(sign(one,fcy(i,j,k)))
                    fym = (one-fracz)*fym + fracz*bY(ii,j,kk)*(x(ii,j,kk)-x(ii,j-1,kk))
                endif 

                fyp = bY(i,j+1,k)*(x(i,j+1,k) - x(i,j,k))
                if (yp) then 
                    fracz = abs(fcy(i,j+1,k))
                    ii = i 
                    kk = k 
                    if(zp) ii = i + int(sign(one,fcy(i,j+1,k)))
                    if(xp) kk = k + int(sign(one,fcy(i,j+1,k)))
                    fyp = (one-fracz)*fyp + fracy*bY(ii,j+1,kk)*(x(ii,j+1,kk)-x(ii,j,kk))
                endif 

                fzm = bZ(i,j,k)*(x(i,j,k) - x(i,j,k-1))
                if (zm) then 
                    fracx = abs(fcz(i,j,k))
                    ii = i
                    jj = j 
                    if(ym) ii = i + int(sign(one,fcz(i,j,k)))
                    if(xm) jj = j + int(sign(one,fcz(i,j,k)))
                    fzm = (one-fracx)*fzm + fracx*bZ(ii,jj,k)*(x(ii,jj,k)-x(ii,jj,k-1))
                endif 

                fzp = bZ(i,j,k+1)*(x(i,j,k+1) - x(i,j,k))
                if (zp) then 
                    fracx = abs(fcz(i+1,j,k))
                    ii = i 
                    jj = j 
                    if(yp) ii = i + int(sign(one,fcz(i,j,k+1)))
                    if(xp) jj = j + int(sign(one,fcz(i,j,k+1)))
                    fzp = (one-fracx)*fzp + fracy*bZ(ii,jj,k+1)*(x(ii,jj,k+1)-x(ii,jj,k))
                endif 


                y(i,j,k) = (one/vfrc(i,j,k))*&
                       (dhx*(apx(i,j,k)*fxm - apx(i+1,j,k)*fxp) + dhy*(apy(i,j,k)*fym & 
                        -apy(i,j+1,k)*fyp) + dhz*(apz(i,j,k)*fzm - apz(i,j,k+1)*fzp))
               if(alpha .ne. zero) then 
                y(i,j,k) = y(i,j,k) + alpha*a(i,j,k)*x(i,j,k)
                call amrex_error("amrex_mlebabeclap_adotx: alpha .ne. 0 todo") 
               endif
              endif 
          enddo
       enddo 
   enddo 
  end subroutine amrex_mlebabeclap_adotx

  subroutine amrex_mlebabeclap_gsrb(lo, hi, phi, hlo, hhi, rhs, rlo, rhi, a, alo, ahi, & 
     bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, & 
     m0, m0lo, m0hi, m2, m2lo, m2hi, & 
     m1, m1lo, m1hi, m3, m3lo, m3hi, & 
     f0, f0lo, f0hi, f2, f2lo, f2hi, & 
     f1, f1lo, f1hi, f3, f3lo, f3hi, & 
     flag, flo, fhi, vfrc, vlo, vhi, & 
     apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
     fcz, czlo, czhi, dxinv, alpha, beta, redblack) & 
     bind(c,name='amrex_mlebabeclap_gsrb') 

 end subroutine amrex_mlebabeclap_gsrb


end module amrex_mlebabeclap_3d_module
