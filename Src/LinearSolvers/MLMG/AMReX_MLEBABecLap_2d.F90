
module amrex_mlebabeclap_2d_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, get_neighbor_cells_int_single
  implicit none

  private
  public :: amrex_mlebabeclap_adotx, amrex_mlebabeclap_gsrb

contains

  subroutine amrex_mlebabeclap_adotx(lo, hi, y, ylo, yhi, x, xlo, xhi, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, flag, flo, fhi, vfrc, vlo, vhi, &
       apx, axlo, axhi, apy, aylo, ayhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
       dxinv, alpha, beta) &
       bind(c,name='amrex_mlebabeclap_adotx')
    integer, dimension(2), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, alo, ahi, bxlo, bxhi, bylo, byhi, &
         flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, cxlo, cxhi, cylo, cyhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), value, intent(in) :: alpha, beta
    real(amrex_real), intent(inout) ::    y( ylo(1): yhi(1), ylo(2): yhi(2))
    real(amrex_real), intent(in   ) ::    x( xlo(1): xhi(1), xlo(2): xhi(2))
    real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2))
    real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2))
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2))
    real(amrex_real), intent(in   ) :: vfrc( vlo(1): vhi(1), vlo(2): vhi(2))
    real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2))
    real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2))
    real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2))

    integer :: i,j
    real(amrex_real) :: dhx, dhy, fxm, fxp, fym, fyp, fracx, fracy

    dhx = beta*dxinv(1)*dxinv(1)
    dhy = beta*dxinv(2)*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (is_covered_cell(flag(i,j))) then
             y(i,j) = 0.d0
          else if (is_regular_cell(flag(i,j))) then
             y(i,j) = alpha*a(i,j)*x(i,j) &
                  - dhx * (bX(i+1,j)*(x(i+1,j) - x(i  ,j))  &
                  &      - bX(i  ,j)*(x(i  ,j) - x(i-1,j))) &
                  - dhy * (bY(i,j+1)*(x(i,j+1) - x(i,j  ))  &
                  &      - bY(i,j  )*(x(i,j  ) - x(i,j-1)))
          else
             fxm = bX(i,j)*(x(i,j)-x(i-1,j))
             if (apx(i,j).ne.0.d0 .and. apx(i,j).ne.1.d0) then
                if (fcx(i,j) .le. 0.d0) then
                   fracy = -fcx(i,j)
                   fxm = (1.d0-fracy)*fxm + fracy*bX(i,j-1)*(x(i,j-1)-x(i-1,j-1))
                else
                   fracy = fcx(i,j)
                   fxm = (1.d0-fracy)*fxm + fracy*bX(i,j+1)*(x(i,j+1)-x(i-1,j+1))
                end if
             end if

             fxp = bX(i+1,j)*(x(i+1,j)-x(i,j))
             if (apx(i+1,j).ne.0.d0 .and. apx(i+1,j).ne.1.d0) then
                if (fcx(i+1,j) .le. 0.d0) then
                   fracy = -fcx(i+1,j)
                   fxp = (1.d0-fracy)*fxp + fracy*bX(i+1,j-1)*(x(i+1,j-1)-x(i,j-1))
                else
                   fracy = fcx(i+1,j)
                   fxp = (1.d0-fracy)*fxp + fracy*bX(i+1,j+1)*(x(i+1,j+1)-x(i,j+1))
                end if
             end if

             fym = bY(i,j)*(x(i,j)-x(i,j-1))
             if (apy(i,j).ne.0.d0 .and. apy(i,j).ne.1.d0) then
                if (fcy(i,j).le.0.d0) then
                   fracx = -fcy(i,j)
                   fym = (1.d0-fracx)*fym + fracx*bY(i-1,j)*(x(i-1,j)-x(i-1,j-1))
                else
                   fracx = fcy(i,j)
                   fym = (1.d0-fracx)*fym + fracx*bY(i+1,j)*(x(i+1,j)-x(i+1,j-1))
                end if
             end if

             fyp = bY(i,j+1)*(x(i,j+1)-x(i,j))
             if (apy(i,j+1).ne.0.d0 .and. apy(i,j+1).ne.1.d0) then
                if (fcy(i,j+1).le.0.d0) then
                   fracx = -fcy(i,j+1)
                   fyp = (1.d0-fracx)*fyp + fracx*bY(i-1,j+1)*(x(i-1,j+1)-x(i-1,j))
                else
                   fracx = fcy(i,j+1)
                   fyp = (1.d0-fracx)*fyp + fracx*bY(i+1,j+1)*(x(i+1,j+1)-x(i+1,j))
                end if
             end if

             y(i,j) = (1.d0/vfrc(i,j)) * &
                  (dhx*(apx(i,j)*fxm-apx(i+1,j)*fxp) + dhy*(apy(i,j)*fym-apy(i,j+1)*fyp))

             if (alpha .ne. 0.d0) then
                y(i,j) = y(i,j) + alpha*a(i,j)*x(i,j)
                call amrex_error("amrex_mlebabeclap_adotx: alpha .ne. 0 todo")
             end if
          end if
       end do
    end do
  end subroutine amrex_mlebabeclap_adotx


  subroutine amrex_mlebabeclap_gsrb(lo, hi, phi, hlo, hhi, rhs, rlo, rhi, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, &
       m0, m0lo, m0hi, m2, m2lo, m2hi, &
       m1, m1lo, m1hi, m3, m3lo, m3hi, &
       f0, f0lo, f0hi, f2, f2lo, f2hi, &
       f1, f1lo, f1hi, f3, f3lo, f3hi, &
       flag, flo, fhi, vfrc, vlo, vhi, &
       apx, axlo, axhi, apy, aylo, ayhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
       dxinv, alpha, beta, redblack) &
       bind(c,name='amrex_mlebabeclap_gsrb')
    integer, dimension(2), intent(in) :: lo, hi, hlo, hhi, rlo, rhi, alo, ahi, bxlo, bxhi, bylo, byhi, &
         m0lo, m0hi, m1lo, m1hi, m2lo, m2hi, m3lo, m3hi, &
         f0lo, f0hi, f1lo, f1hi, f2lo, f2hi, f3lo, f3hi, &
         flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, cxlo, cxhi, cylo, cyhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), value, intent(in) :: alpha, beta
    integer, value, intent(in) :: redblack
    real(amrex_real), intent(inout) ::  phi( hlo(1): hhi(1), hlo(2): hhi(2))
    real(amrex_real), intent(in   ) ::  rhs( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) ::    a( alo(1): ahi(1), alo(2): ahi(2))
    real(amrex_real), intent(in   ) ::   bx(bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(amrex_real), intent(in   ) ::   by(bylo(1):byhi(1),bylo(2):byhi(2))
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

    integer :: i,j,ioff
    real(amrex_real) :: cf0, cf1, cf2, cf3, delta, gamma, rho
    real(amrex_real) :: dhx, dhy, fxm, fxp, fym, fyp, fracx, fracy
    real(amrex_real) :: sxm, sxp, sym, syp

    dhx = beta*dxinv(1)*dxinv(1)
    dhy = beta*dxinv(2)*dxinv(2)

    if (alpha.ne.0.d0) call amrex_error("amrex_mlebabeclap_gsrb: todo")

    do j = lo(2), hi(2)
       ioff = mod(lo(1)+j+redblack,2)
       do i = lo(1)+ioff, hi(1), 2

          if (is_covered_cell(flag(i,j))) then
             phi(i,j) = 0.d0
          else
               cf0 = merge(f0(lo(1),j), 0.0D0, &
                    (i .eq. lo(1)) .and. (m0(lo(1)-1,j).gt.0))
               cf1 = merge(f1(i,lo(2)), 0.0D0, &
                    (j .eq. lo(2)) .and. (m1(i,lo(2)-1).gt.0))
               cf2 = merge(f2(hi(1),j), 0.0D0, &
                    (i .eq. hi(1)) .and. (m2(hi(1)+1,j).gt.0))
               cf3 = merge(f3(i,hi(2)), 0.0D0, &
                    (j .eq. hi(2)) .and. (m3(i,hi(2)+1).gt.0))

               delta = dhx*(bX(i,j)*cf0 + bX(i+1,j)*cf2) &
                    +  dhy*(bY(i,j)*cf1 + bY(i,j+1)*cf3)

               if (is_regular_cell(flag(i,j))) then

                  gamma = alpha*a(i,j) &
                       + dhx * (bX(i+1,j) + bX(i,j)) &
                       + dhy * (bY(i,j+1) + bY(i,j))

                  rho =  dhx * (bX(i+1,j)*phi(i+1,j) + bX(i,j)*phi(i-1,j)) &
                       + dhy * (bY(i,j+1)*phi(i,j+1) + bY(i,j)*phi(i,j-1))

               else

                  fxm = -bX(i,j)*phi(i-1,j)
                  sxm =  bX(i,j)
                  if (apx(i,j).ne.0.d0 .and. apx(i,j).ne.1.d0) then
                     if (fcx(i,j) .le. 0.d0) then
                        fracy = -fcx(i,j)
                        fxm = (1.d0-fracy)*fxm + fracy*bX(i,j-1)*(phi(i,j-1)-phi(i-1,j-1))
                        sxm = (1.d0-fracy)*sxm
                     else
                        fracy = fcx(i,j)
                        fxm = (1.d0-fracy)*fxm + fracy*bX(i,j+1)*(phi(i,j+1)-phi(i-1,j+1))
                        sxm = (1.d0-fracy)*sxm
                     end if
                  end if
                  
                  fxp =  bX(i+1,j)*phi(i+1,j)
                  sxp = -bX(i+1,j)
                  if (apx(i+1,j).ne.0.d0 .and. apx(i+1,j).ne.1.d0) then
                     if (fcx(i+1,j) .le. 0.d0) then
                        fracy = -fcx(i+1,j)
                        fxp = (1.d0-fracy)*fxp + fracy*bX(i+1,j-1)*(phi(i+1,j-1)-phi(i,j-1))
                        sxp = (1.d0-fracy)*sxp
                     else
                        fracy = fcx(i+1,j)
                        fxp = (1.d0-fracy)*fxp + fracy*bX(i+1,j+1)*(phi(i+1,j+1)-phi(i,j+1))
                        sxp = (1.d0-fracy)*sxp
                     end if
                  end if
                  
                  fym = -bY(i,j)*phi(i,j-1)
                  sym =  bY(i,j)
                  if (apy(i,j).ne.0.d0 .and. apy(i,j).ne.1.d0) then
                     if (fcy(i,j).le.0.d0) then
                        fracx = -fcy(i,j)
                        fym = (1.d0-fracx)*fym + fracx*bY(i-1,j)*(phi(i-1,j)-phi(i-1,j-1))
                        sym = (1.d0-fracx)*sym
                     else
                        fracx = fcy(i,j)
                        fym = (1.d0-fracx)*fym + fracx*bY(i+1,j)*(phi(i+1,j)-phi(i+1,j-1))
                        sym = (1.d0-fracx)*sym
                     end if
                  end if
                  
                  fyp =  bY(i,j+1)*phi(i,j+1)
                  syp = -bY(i,j+1)
                  if (apy(i,j+1).ne.0.d0 .and. apy(i,j+1).ne.1.d0) then
                     if (fcy(i,j+1).le.0.d0) then
                        fracx = -fcy(i,j+1)
                        fyp = (1.d0-fracx)*fyp + fracx*bY(i-1,j+1)*(phi(i-1,j+1)-phi(i-1,j))
                        syp = (1.d0-fracx)*syp
                     else
                        fracx = fcy(i,j+1)
                        fyp = (1.d0-fracx)*fyp + fracx*bY(i+1,j+1)*(phi(i+1,j+1)-phi(i+1,j))
                        syp = (1.d0-fracx)*syp
                     end if
                  end if

                  gamma = alpha*a(i,j) + (1.d0/vfrc(i,j)) * &
                       (dhx*(apx(i,j)*sxm-apx(i+1,j)*sxp) + dhy*(apy(i,j)*sym-apy(i,j+1)*syp))

                  rho = -(1.d0/vfrc(i,j)) * &
                       (dhx*(apx(i,j)*fxm-apx(i+1,j)*fxp) + dhy*(apy(i,j)*fym-apy(i,j+1)*fyp))
               end if

               phi(i,j) = (rhs(i,j) + rho - phi(i,j)*delta) / (gamma - delta)
            end if
       end do
    end do

  end subroutine amrex_mlebabeclap_gsrb

end module amrex_mlebabeclap_2d_module
