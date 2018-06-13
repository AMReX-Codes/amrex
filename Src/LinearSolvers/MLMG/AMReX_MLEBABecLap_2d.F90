
module amrex_mlebabeclap_2d_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, get_neighbor_cells_int_single
  implicit none

  private
  public :: amrex_mlebabeclap_adotx

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
                call amrex_error("amrex_mlebabeclap_adotx: alpha .ne. 0 todo")
             end if
          end if
       end do
    end do
  end subroutine amrex_mlebabeclap_adotx

end module amrex_mlebabeclap_2d_module
