
module amrex_eb2_2d_moudle

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module, only : zero, one, half, sixth, eighth
  implicit none

  integer, parameter :: regular = 0
  integer, parameter :: covered = 1
  integer, parameter :: irregular = 2
  integer, parameter :: unknown = 3

  real(amrex_real), private, parameter :: small = 1.d-14

  private
  public :: amrex_eb2_gfab_build_types, amrex_eb2_build_faces, amrex_eb2_build_cells

contains

  subroutine amrex_eb2_gfab_build_types (lo, hi, s, slo, shi, cell, clo, chi, &
       fx, fxlo, fxhi, fy, fylo, fyhi) bind(c,name='amrex_eb2_gfab_build_types')
    integer, dimension(2), intent(in) :: lo, hi, slo, shi, clo, chi, fxlo, fxhi, fylo, fyhi
    real(amrex_real), intent(in   ) ::    s( slo(1): shi(1), slo(2): shi(2))
    integer(c_int)  , intent(inout) :: cell( clo(1): chi(1), clo(2): chi(2))
    integer(c_int)  , intent(inout) ::   fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    integer(c_int)  , intent(inout) ::   fy(fylo(1):fyhi(1),fylo(2):fyhi(2))

    integer :: i,j

    stop "todo"

    do    j = lo(2)-2, hi(2)+2
       do i = lo(1)-2, hi(1)+2
          if (       s(i,j  ).ge.zero .and. s(i+1,j  ).ge.zero &
               .and. s(i,j+1).ge.zero .and. s(i+1,j+1).ge.zero) then
             cell(i,j) = covered
          else if (  s(i,j  ).lt.zero .and. s(i+1,j  ).lt.zero &
               .and. s(i,j+1).lt.zero .and. s(i+1,j+1).lt.zero) then
             cell(i,j) = regular
          else
             cell(i,j) = irregular
          end if
       end do
    end do

    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+2
          if (s(i,j).ge.zero .and. s(i,j+1).ge.zero) then
             fx(i,j) = covered
          else if (s(i,j).lt.zero .and. s(i,j+1).lt.zero) then
             fx(i,j) = regular
          else
             fx(i,j) = irregular
          end if
       end do
    end do

    do    j = lo(2)-1, hi(2)+2
       do i = lo(1)-1, hi(1)+1
          if (s(i,j).ge.zero .and. s(i+1,j).ge.zero) then
             fy(i,j) = covered
          else if (s(i,j).lt.zero .and. s(i+1,j).lt.zero) then
             fy(i,j) = regular
          else
             fy(i,j) = irregular
          end if
       end do
    end do

  end subroutine amrex_eb2_gfab_build_types


  subroutine amrex_eb2_build_faces (lo, hi, cell, clo, chi, &
       fx, fxlo, fxhi, fy, fylo, fyhi, levset, slo, shi,&
       interx, ixlo, ixhi, intery, iylo, iyhi, &
       apx, axlo, axhi, apy, aylo, ayhi, &
       fcx, cxlo, cxhi, fcy, cylo, cyhi, &
       dx, dxinv, problo) &
       bind(c,name='amrex_eb2_build_faces')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, slo, shi, &
         fxlo, fxhi, fylo, fyhi, ixlo, ixhi, iylo, iyhi, axlo, axhi, aylo, ayhi, &
         cxlo, cxhi, cylo, cyhi
    real(amrex_real), dimension(2) :: dx, dxinv, problo
    integer(c_int)  , intent(inout) ::   cell( clo(1): chi(1), clo(2): chi(2))
    integer(c_int)  , intent(inout) ::     fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    integer(c_int)  , intent(inout) ::     fy(fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(amrex_real), intent(in   ) :: levset( slo(1): shi(1), slo(2): shi(2))
    real(amrex_real), intent(in   ) :: interx(ixlo(1):ixhi(1),ixlo(2):ixhi(2))
    real(amrex_real), intent(in   ) :: intery(iylo(1):iyhi(1),iylo(2):iyhi(2))
    real(amrex_real), intent(inout) ::    apx(axlo(1):axhi(1),axlo(2):axhi(2))
    real(amrex_real), intent(inout) ::    apy(aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(amrex_real), intent(inout) ::    fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),2)
    real(amrex_real), intent(inout) ::    fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),2)

    integer :: i,j, ncuts

    ! x-face
    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+2
          if (fx(i,j) .eq. regular) then
             apx(i,j) = one
             fcx(i,j,:) = zero
          else if (fx(i,j) .eq. covered) then
             apx(i,j) = zero
             fcx(i,j,:) = zero
          else
             if (levset(i,j) .lt. zero) then
                apx(i,j) = (interx(i,j)-(problo(2)+j*dx(2)))*dxinv(2)
                fcx(i,j,1) = zero
                fcx(i,j,2) = half*apx(i,j) - half
             else
                apx(i,j) = one - (interx(i,j)-(problo(2)+j*dx(2)))*dxinv(2)
                fcx(i,j,1) = zero
                fcx(i,j,2) = half - half*apx(i,j)
             end if

             if (apx(i,j) .gt. one-small) then
                apx(i,j) = one
                fcx(i,j,:) = zero
                fx(i,j) = regular
             else if (apx(i,j) .lt. small) then
                apx(i,j) = zero
                fcx(i,j,:) = zero
                fx(i,j) = covered
             end if
          end if
       end do
    end do

    ! y-face
    do    j = lo(2)-1, hi(2)+2
       do i = lo(1)-1, hi(1)+1
          if (fy(i,j) .eq. regular) then
             apy(i,j) = one
             fcy(i,j,:) = zero
          else if (fy(i,j) .eq. covered) then
             apy(i,j) = zero
             fcy(i,j,:) = zero
          else
             if (levset(i,j) .lt. zero) then
                apy(i,j) = (intery(i,j)-(problo(1)+i*dx(1)))*dxinv(1)
                fcy(i,j,1) = half*apy(i,j) - half
                fcy(i,j,2) = zero
             else
                apy(i,j) = one - (intery(i,j)-(problo(1)+i*dx(1)))*dxinv(1)
                fcy(i,j,1) = half - half*apy(i,j)
                fcy(i,j,2) = zero
             end if

             if (apy(i,j) .gt. one-small) then
                apy(i,j) = one
                fcy(i,j,:) = zero
                fy(i,j) = regular
             else if (fy(i,j) .eq. covered) then
                apy(i,j) = zero
                fcy(i,j,:) = zero
                fy(i,j) = covered
             end if
          end if
       end do
    end do

    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+1
          if (cell(i,j) .eq. irregular) then
             if ( fx(i,j).eq.regular .and. fx(i+1,j).eq.regular .and. &
                  fy(i,j).eq.regular .and. fy(i,j+1).eq.regular ) then
                cell(i,j) = regular
             else if ( fx(i,j).eq.covered .and. fx(i+1,j).eq.covered .and. &
                  &    fy(i,j).eq.covered .and. fy(i,j+1).eq.covered ) then
                cell(i,j) = covered
             else
                ncuts = 0
                if (fx(i,j).eq.irregular) ncuts = ncuts+1
                if (fx(i+1,j).eq.irregular) ncuts = ncuts+1
                if (fy(i,j).eq.irregular) ncuts = ncuts+1
                if (fy(i,j+1).eq.irregular) ncuts = ncuts+1
                if (ncuts .ne. 2) then
                   call amrex_error("amrex_eb2_build_faces: wrong number of cuts")
                end if
             end if
          end if
       end do
    end do

  end subroutine amrex_eb2_build_faces


  subroutine amrex_eb2_build_cells (lo, hi, cell, clo, chi, &
       fx, fxlo, fxhi, fy, fylo, fyhi, &
       apx, axlo, axhi, apy, aylo, ayhi, &
       vfrac, vlo, vhi, vcent, tlo, thi, barea, alo, ahi, &
       bcent, blo, bhi, bnorm, mlo, mhi) &
       bind(c,name='amrex_eb2_build_cells')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, &
         fxlo, fxhi, fylo, fyhi, axlo, axhi, aylo, ayhi, &
         vlo, vhi, tlo, thi, alo, ahi, blo, bhi, mlo, mhi
    integer(c_int)  , intent(inout) ::  cell( clo(1): chi(1), clo(2): chi(2))
    integer(c_int)  , intent(inout) ::    fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    integer(c_int)  , intent(inout) ::    fy(fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(amrex_real), intent(inout) ::   apx(axlo(1):axhi(1),axlo(2):axhi(2))
    real(amrex_real), intent(inout) ::   apy(aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(amrex_real), intent(inout) :: vfrac( vlo(1): vhi(1), vlo(2): vhi(2))
    real(amrex_real), intent(inout) :: vcent( tlo(1): thi(1), tlo(2): thi(2),2)
    real(amrex_real), intent(inout) :: barea( alo(1): ahi(1), alo(2): ahi(2))
    real(amrex_real), intent(inout) :: bcent( blo(1): bhi(1), blo(2): bhi(2),2)
    real(amrex_real), intent(inout) :: bnorm( mlo(1): mhi(1), mlo(2): mhi(2),2)

    integer :: i,j
    real(amrex_real) :: axm, axp, aym, ayp, apnorm, apnorminv
    real(amrex_real) :: nx, ny, x_ym, x_yp, y_xm, y_xp, aa, af1, af2
    real(amrex_real) :: dx, dx2, dx3, dx4, dy, dy2, dy3, dy4, nxabs, nyabs, signx, signy
    real(amrex_real), parameter :: tiny = 1.d-15
    real(amrex_real), parameter :: almostone = 1.d0-1.d-15

    ! xxxx clean up local variables

    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+1
          if (cell(i,j) .eq. regular) then
             vfrac(i,j) = one
             vcent(i,j,:) = zero
             barea(i,j) = zero
             bcent(i,j,:) = zero
             bnorm(i,j,:) = zero
          else if (cell(i,j) .eq. covered) then
             vfrac(i,j) = zero
             vcent(i,j,:) = zero
             barea(i,j) = zero
             bcent(i,j,:) = zero
             bnorm(i,j,:) = zero
          else
             axm = apx(i,j)
             axp = apx(i+1,j)
             aym = apy(i,j)
             ayp = apy(i,j+1)
             apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2)
             apnorminv = one/apnorm
             nx = (axm-axp) * apnorminv  ! pointing to the wall
             ny = (aym-ayp) * apnorminv

             nxabs = abs(nx)
             nyabs = abs(ny)

             if (nx .gt. zero) then
                x_ym = -half + aym
                x_yp = -half + ayp
                signx = one
             else
                x_ym = half - aym
                x_yp = half - ayp
                signx = -one
             end if

             if (ny .gt. zero) then
                y_xm = -half + axm
                y_xp = -half + axp
                signy = one
             else
                y_xm = half - axm
                y_xp = half - axp
                signy = -one
             end if

             barea(i,j) = nx*(axm-axp) + ny*(aym-ayp)
             bcent(i,j,1) = half*(x_ym+x_yp)
             bcent(i,j,2) = half*(y_xm+y_xp)
             bnorm(i,j,1) = nx
             bnorm(i,j,2) = ny

             if (nxabs .lt. tiny .or. nyabs .gt. almostone) then
                barea(i,j) = one
                bcent(i,j,1) = zero
                bnorm(i,j,1) = zero
                bnorm(i,j,2) = signy
                vfrac(i,j) = half*(axm+axp)
                vcent(i,j,1) = zero
                vcent(i,j,2) = (eighth*(ayp-aym) + ny*half*bcent(i,j,2)**2)/vfrac(i,j)
             else if (nyabs .lt. tiny .or. nxabs .gt. almostone) then
                barea(i,j) = one
                bcent(i,j,2) = zero
                bnorm(i,j,1) = signx
                bnorm(i,j,2) = zero
                vfrac(i,j) = half*(aym+ayp)
                vcent(i,j,1) = (eighth*(axp-axm) + nx*half*bcent(i,j,1)**2)/vfrac(i,j)
                vcent(i,j,2) = zero
             else
                aa = nxabs/ny
                dx  = x_ym - x_yp
                dx2 = dx * (x_ym + x_yp)
                dx3 = dx * (x_ym**2 + x_ym*x_yp + x_yp**2)
                dx4 = dx * (x_ym + x_yp) * (x_ym**2 + x_yp**2)
                af1 = half*(axm+axp) + aa*half*dx2
                vcent(i,j,1) = eighth*(axp-axm) + aa*sixth*dx3

                aa = nyabs/nx
                dy = y_xm - y_xp
                dy2 = dy * (y_xm + y_xp)
                dy3 = dy * (y_xm**2 + y_xm*y_xp + y_xp**2)
                dy4 = dy * (y_xm + y_xp) * (y_xm**2 + y_xp**2)
                af2 = half*(aym+ayp) + aa*half*dy2
                vcent(i,j,2) = eighth*(ayp-aym) + aa*sixth*dy3

                vfrac(i,j) = half*(af1+af2)
                if (vfrac(i,j) .gt. one-small) then
                   vfrac(i,j) = one
                   vcent(i,j,1) = zero
                   vcent(i,j,2) = zero
                else if (vfrac(i,j) .lt. small) then
                   vfrac(i,j) = zero
                   vcent(i,j,1) = zero
                   vcent(i,j,2) = zero
                else
                   vcent(i,j,1) = vcent(i,j,1)*(one/vfrac(i,j))
                   vcent(i,j,2) = vcent(i,j,2)*(one/vfrac(i,j))
                   vcent(i,j,1) = min(max(vcent(i,j,1),-half),half)
                   vcent(i,j,2) = min(max(vcent(i,j,2),-half),half)
                end if
             end if
          end if
       end do
    end do
  end subroutine amrex_eb2_build_cells

end module amrex_eb2_2d_moudle
