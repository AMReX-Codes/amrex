
module amrex_eb2_2d_moudle

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module, only : zero, one, half, fourth, sixth, eighth
  use amrex_ebcellflag_module, only : regular_cell => regular, covered_cell => covered, &
       is_regular_cell, is_single_valued_cell, is_covered_cell, get_cell_type, get_neighbor_cells, &
       set_regular_cell, set_single_valued_cell, set_covered_cell, &
       set_neighbor, clear_neighbor, is_neighbor
  implicit none

  integer, parameter :: regular = 0
  integer, parameter :: covered = 1
  integer, parameter :: irregular = 2
  integer, parameter :: unknown = 3

  real(amrex_real), private, parameter :: small = 1.d-14

  private
  public :: amrex_eb2_gfab_build_types, amrex_eb2_build_faces, amrex_eb2_build_cells, &
       amrex_eb2_coarsen_from_fine, amrex_eb2_build_cellflag_from_ap, amrex_eb2_check_mvmc

contains

  subroutine amrex_eb2_gfab_build_types (lo, hi, s, slo, shi, cell, clo, chi, &
       fx, fxlo, fxhi, fy, fylo, fyhi) bind(c,name='amrex_eb2_gfab_build_types')
    integer, dimension(2), intent(in) :: lo, hi, slo, shi, clo, chi, fxlo, fxhi, fylo, fyhi
    real(amrex_real), intent(in   ) ::    s( slo(1): shi(1), slo(2): shi(2))
    integer(c_int)  , intent(inout) :: cell( clo(1): chi(1), clo(2): chi(2))
    integer(c_int)  , intent(inout) ::   fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    integer(c_int)  , intent(inout) ::   fy(fylo(1):fyhi(1),fylo(2):fyhi(2))

    integer :: i,j

    do    j = lo(2)-2, hi(2)+2
       do i = lo(1)-2, hi(1)+2
          if (       s(i,j  ).ge.zero .and. s(i+1,j  ).ge.zero &
               .and. s(i,j+1).ge.zero .and. s(i+1,j+1).ge.zero) then
             call set_covered_cell(cell(i,j))
          else if (  s(i,j  ).lt.zero .and. s(i+1,j  ).lt.zero &
               .and. s(i,j+1).lt.zero .and. s(i+1,j+1).lt.zero) then
             call set_regular_cell(cell(i,j))
          else
             call set_single_valued_cell(cell(i,j))
          end if
       end do
    end do

    do    j = lo(2)-2, hi(2)+2
       do i = lo(1)-2, hi(1)+3
          if (s(i,j).ge.zero .and. s(i,j+1).ge.zero) then
             fx(i,j) = covered
          else if (s(i,j).lt.zero .and. s(i,j+1).lt.zero) then
             fx(i,j) = regular
          else
             fx(i,j) = irregular
          end if
       end do
    end do

    do    j = lo(2)-2, hi(2)+3
       do i = lo(1)-2, hi(1)+2
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
    real(amrex_real), intent(inout) ::    fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2))
    real(amrex_real), intent(inout) ::    fcy(cylo(1):cyhi(1),cylo(2):cyhi(2))

    integer :: i,j, ncuts

    ! x-face
    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+2
          if (fx(i,j) .eq. regular) then
!             apx(i,j) = one
!             fcx(i,j) = zero
          else if (fx(i,j) .eq. covered) then
             apx(i,j) = zero
             fcx(i,j) = zero
          else
             if (levset(i,j) .lt. zero) then
                apx(i,j) = (interx(i,j)-(problo(2)+j*dx(2)))*dxinv(2)
                fcx(i,j) = half*apx(i,j) - half
             else
                apx(i,j) = one - (interx(i,j)-(problo(2)+j*dx(2)))*dxinv(2)
                fcx(i,j) = half - half*apx(i,j)
             end if

             if (apx(i,j) .gt. one-small) then
                apx(i,j) = one
                fcx(i,j) = zero
                fx(i,j) = regular
             else if (apx(i,j) .lt. small) then
                apx(i,j) = zero
                fcx(i,j) = zero
                fx(i,j) = covered
             end if
          end if
       end do
    end do

    ! y-face
    do    j = lo(2)-1, hi(2)+2
       do i = lo(1)-1, hi(1)+1
          if (fy(i,j) .eq. regular) then
!             apy(i,j) = one
!             fcy(i,j) = zero
          else if (fy(i,j) .eq. covered) then
             apy(i,j) = zero
             fcy(i,j) = zero
          else
             if (levset(i,j) .lt. zero) then
                apy(i,j) = (intery(i,j)-(problo(1)+i*dx(1)))*dxinv(1)
                fcy(i,j) = half*apy(i,j) - half
             else
                apy(i,j) = one - (intery(i,j)-(problo(1)+i*dx(1)))*dxinv(1)
                fcy(i,j) = half - half*apy(i,j)
             end if

             if (apy(i,j) .gt. one-small) then
                apy(i,j) = one
                fcy(i,j) = zero
                fy(i,j) = regular
             else if (apy(i,j) .lt. small) then
                apy(i,j) = zero
                fcy(i,j) = zero
                fy(i,j) = covered
             end if
          end if
       end do
    end do

    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+1
          if (is_single_valued_cell(cell(i,j))) then
             if ( fx(i,j).eq.regular .and. fx(i+1,j).eq.regular .and. &
                  fy(i,j).eq.regular .and. fy(i,j+1).eq.regular ) then
                call set_regular_cell(cell(i,j))
             else if ( fx(i,j).eq.covered .and. fx(i+1,j).eq.covered .and. &
                  &    fy(i,j).eq.covered .and. fy(i,j+1).eq.covered ) then
                call set_covered_cell(cell(i,j))
             else
                ncuts = 0
                if (fx(i,j).eq.irregular) ncuts = ncuts+1
                if (fx(i+1,j).eq.irregular) ncuts = ncuts+1
                if (fy(i,j).eq.irregular) ncuts = ncuts+1
                if (fy(i,j+1).eq.irregular) ncuts = ncuts+1
                if (ncuts .gt. 2) then
                   call amrex_error("amrex_eb2_build_faces: more than 2 cuts not supported")
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

    integer :: i,j, flg
    real(amrex_real) :: axm, axp, aym, ayp, apnorm, apnorminv
    real(amrex_real) :: nx, ny, x_ym, x_yp, y_xm, y_xp, aa, af1, af2
    real(amrex_real) :: dx, dx2, dx3, dx4, dy, dy2, dy3, dy4, nxabs, nyabs, signx, signy
    real(amrex_real), parameter :: tiny = 1.d-15
    real(amrex_real), parameter :: almostone = 1.d0-1.d-15

    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+1
          if (is_regular_cell(cell(i,j))) then
!             vfrac(i,j) = one
!             vcent(i,j,:) = zero
!             barea(i,j) = zero
!             bcent(i,j,:) = -one
!             bnorm(i,j,:) = zero
          else if (is_covered_cell(cell(i,j))) then
             vfrac(i,j) = zero
!             vcent(i,j,:) = zero
!             barea(i,j) = zero
!             bcent(i,j,:) = -one
!             bnorm(i,j,:) = zero
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

             ! remove small cells
             if (vfrac(i,j) < small) then
                vfrac(i,j) = zero
                vcent(i,j,:) = zero
                barea(i,j) = zero
                bcent(i,j,:) = -one
                bnorm(i,j,:) = zero
                call set_covered_cell(cell(i,j))
             end if
          end if
       end do
    end do

    ! fix faces for small cells
    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)+1
          if (vfrac(i-1,j) < small .or. vfrac(i,j) < small) then
             fx(i,j) = covered
             apx(i,j) = zero
          end if
       end do
    end do
    !
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-1, hi(1)+1
          if (vfrac(i,j-1) < small .or. vfrac(i,j) < small) then
             fy(i,j) = covered
             apy(i,j) = zero
          end if
       end do
    end do

    ! Build neighbors.  By default, all neighbors are already set.
    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+1
          flg = cell(i,j)

          if (fx(i  ,j  ).eq.covered) flg = clear_neighbor(flg, -1,  0)
          if (fx(i+1,j  ).eq.covered) flg = clear_neighbor(flg,  1,  0)
          if (fy(i  ,j  ).eq.covered) flg = clear_neighbor(flg,  0, -1)
          if (fy(i  ,j+1).eq.covered) flg = clear_neighbor(flg,  0,  1)

          if (fx(i,j).ne.covered .and. fy(i-1,j).ne.covered) then
          else if (fx(i,j-1).ne.covered .and. fy(i,j).ne.covered) then
          else
             flg = clear_neighbor(flg,-1,-1)
          end if

          if (fx(i+1,j).ne.covered .and. fy(i+1,j).ne.covered) then
          else if (fx(i+1,j-1).ne.covered .and. fy(i,j).ne.covered) then
          else
             flg = clear_neighbor(flg,1,-1)
          end if

          if (fx(i,j).ne.covered .and. fy(i-1,j+1).ne.covered) then
          else if (fx(i,j+1).ne.covered .and. fy(i,j+1).ne.covered) then
          else
             flg = clear_neighbor(flg,-1,1)
          end if

          if (fx(i+1,j).ne.covered .and. fy(i+1,j+1).ne.covered) then
          else if (fx(i+1,j+1).ne.covered .and. fy(i,j+1).ne.covered) then
          else
             flg = clear_neighbor(flg,1,1)
          end if

          cell(i,j) = flg
       end do
    end do

  end subroutine amrex_eb2_build_cells


  subroutine amrex_eb2_coarsen_from_fine (lo, hi, xlo, xhi, ylo, yhi, &
       cvol, cvlo, cvhi, fvol, fvlo, fvhi, ccent, cclo, cchi, fcent, fclo, fchi, &
       cba, cbalo, cbahi, fba, fbalo, fbahi, cbc, cbclo, cbchi, fbc, fbclo, fbchi, &
       cbn, cbnlo, cbnhi, fbn, fbnlo, fbnhi, capx, caxlo, caxhi, fapx, faxlo, faxhi, &
       capy, caylo, cayhi, fapy, faylo, fayhi, &
       cfcx, cfxlo, cfxhi, ffcx, ffxlo, ffxhi, cfcy, cfylo, cfyhi, ffcy, ffylo, ffyhi, &
       cflag, cflo, cfhi, fflag, fflo, ffhi, ierr) &
       bind(c, name='amrex_eb2_coarsen_from_fine')
    integer, dimension(2), intent(in) :: lo, hi, xlo, xhi, ylo, yhi, &
         cvlo, cvhi,  fvlo, fvhi, cclo, cchi, fclo, fchi, &
         cbalo, cbahi, fbalo, fbahi, cbclo, cbchi, fbclo, fbchi, &
         cbnlo, cbnhi, fbnlo, fbnhi, caxlo, caxhi, faxlo, faxhi, &
         caylo, cayhi, faylo, fayhi, &
         cfxlo, cfxhi, ffxlo, ffxhi, cfylo, cfyhi, ffylo, ffyhi, &
         cflo, cfhi, fflo, ffhi
    integer, intent(inout) :: ierr
    real(amrex_real), intent(inout) :: cvol ( cvlo(1): cvhi(1), cvlo(2): cvhi(2))
    real(amrex_real), intent(in   ) :: fvol ( fvlo(1): fvhi(1), fvlo(2): fvhi(2))
    real(amrex_real), intent(inout) :: ccent( cclo(1): cchi(1), cclo(2): cchi(2),2)
    real(amrex_real), intent(in   ) :: fcent( fclo(1): fchi(1), fclo(2): fchi(2),2)
    real(amrex_real), intent(inout) :: cba  (cbalo(1):cbahi(1),cbalo(2):cbahi(2))
    real(amrex_real), intent(in   ) :: fba  (fbalo(1):fbahi(1),fbalo(2):fbahi(2))
    real(amrex_real), intent(inout) :: cbc  (cbclo(1):cbchi(1),cbclo(2):cbchi(2),2)
    real(amrex_real), intent(in   ) :: fbc  (fbclo(1):fbchi(1),fbclo(2):fbchi(2),2)
    real(amrex_real), intent(inout) :: cbn  (cbnlo(1):cbnhi(1),cbnlo(2):cbnhi(2),2)
    real(amrex_real), intent(in   ) :: fbn  (fbnlo(1):fbnhi(1),fbnlo(2):fbnhi(2),2)
    real(amrex_real), intent(inout) :: capx (caxlo(1):caxhi(1),caxlo(2):caxhi(2))
    real(amrex_real), intent(in   ) :: fapx (faxlo(1):faxhi(1),faxlo(2):faxhi(2))
    real(amrex_real), intent(inout) :: capy (caylo(1):cayhi(1),caylo(2):cayhi(2))
    real(amrex_real), intent(in   ) :: fapy (faylo(1):fayhi(1),faylo(2):fayhi(2))
    real(amrex_real), intent(inout) :: cfcx (cfxlo(1):cfxhi(1),cfxlo(2):cfxhi(2))
    real(amrex_real), intent(in   ) :: ffcx (ffxlo(1):ffxhi(1),ffxlo(2):ffxhi(2))
    real(amrex_real), intent(inout) :: cfcy (cfylo(1):cfyhi(1),cfylo(2):cfyhi(2))
    real(amrex_real), intent(in   ) :: ffcy (ffylo(1):ffyhi(1),ffylo(2):ffyhi(2))
    integer         , intent(inout) :: cflag( cflo(1): cfhi(1), cflo(2): cfhi(2))
    integer         , intent(in   ) :: fflag( fflo(1): ffhi(1), fflo(2): ffhi(2))

    integer :: i,j, ii,jj, ftype(2,2)
    real(amrex_real) :: cvolinv, cbainv, nx, ny, nfac, apinv

    do    j = lo(2), hi(2)
       jj = j*2
       do i = lo(1), hi(1)
          ii = i*2

          ftype = get_cell_type(fflag(ii:ii+1,jj:jj+1))
          if (all(ftype.eq.regular_cell)) then
             ! nothing to do
          else if (all(ftype.eq.covered_cell)) then
             call set_covered_cell(cflag(i,j))
             cvol(i,j) = zero
          else
             call set_single_valued_cell(cflag(i,j))

             cvol(i,j) = fourth*sum(fvol(ii:ii+1,jj:jj+1))
             cvolinv = one/cvol(i,j)
                
             ccent(i,j,1) = fourth*cvolinv* &
                  ( fvol(ii  ,jj  )*(half*fcent(ii  ,jj  ,1)-fourth) &
                  + fvol(ii+1,jj  )*(half*fcent(ii+1,jj  ,1)+fourth) &
                  + fvol(ii  ,jj+1)*(half*fcent(ii  ,jj+1,1)-fourth) &
                  + fvol(ii+1,jj+1)*(half*fcent(ii+1,jj+1,1)+fourth) )
             ccent(i,j,2) = fourth*cvolinv* &
                  ( fvol(ii  ,jj  )*(half*fcent(ii  ,jj  ,2)-fourth) &
                  + fvol(ii+1,jj  )*(half*fcent(ii+1,jj  ,2)-fourth) &
                  + fvol(ii  ,jj+1)*(half*fcent(ii  ,jj+1,2)+fourth) &
                  + fvol(ii+1,jj+1)*(half*fcent(ii+1,jj+1,2)+fourth) )
                
             cba(i,j) = half*sum(fba(ii:ii+1,jj:jj+1))
             cbainv = one/cba(i,j)
                
             cbc(i,j,1) = half*cbainv* &
                  ( fba(ii  ,jj  )*(half*fbc(ii  ,jj  ,1)-fourth) &
                  + fba(ii+1,jj  )*(half*fbc(ii+1,jj  ,1)+fourth) &
                  + fba(ii  ,jj+1)*(half*fbc(ii  ,jj+1,1)-fourth) &
                  + fba(ii+1,jj+1)*(half*fbc(ii+1,jj+1,1)+fourth) )
             cbc(i,j,2) = half*cbainv* &
                  ( fba(ii  ,jj  )*(half*fbc(ii  ,jj  ,2)-fourth) &
                  + fba(ii+1,jj  )*(half*fbc(ii+1,jj  ,2)-fourth) &
                  + fba(ii  ,jj+1)*(half*fbc(ii  ,jj+1,2)+fourth) &
                  + fba(ii+1,jj+1)*(half*fbc(ii+1,jj+1,2)+fourth) )
                
             nx =   fbn(ii  ,jj  ,1)*fba(ii  ,jj  ) &
                  + fbn(ii+1,jj  ,1)*fba(ii+1,jj  ) &
                  + fbn(ii  ,jj+1,1)*fba(ii  ,jj+1) &
                  + fbn(ii+1,jj+1,1)*fba(ii+1,jj+1)
             ny =   fbn(ii  ,jj  ,2)*fba(ii  ,jj  ) &
                  + fbn(ii+1,jj  ,2)*fba(ii+1,jj  ) &
                  + fbn(ii  ,jj+1,2)*fba(ii  ,jj+1) &
                  + fbn(ii+1,jj+1,2)*fba(ii+1,jj+1)
             if (nx.eq.zero .and. ny.eq.zero) then
                ierr = 1
                return
             end if
             nfac = one/sqrt(nx*nx+ny*ny)
             cbn(i,j,1) = nx*nfac
             cbn(i,j,2) = ny*nfac
          end if
       end do
    end do

    do    j = xlo(2), xhi(2)
       jj = j*2
       do i = xlo(1), xhi(1)
          ii = i*2
          
          capx(i,j) = half*(fapx(ii,jj)+fapx(ii,jj+1))
          if (capx(i,j) .ne. zero) then
             apinv = one/capx(i,j)
             cfcx(i,j) = half*apinv* &
                  ( fapx(ii,jj  )*(half*ffcx(ii,jj  )-fourth) &
                  + fapx(ii,jj+1)*(half*ffcx(ii,jj+1)+fourth) )
          end if
       end do
    end do

    do    j = ylo(2), yhi(2)
       jj = j*2
       do i = ylo(1), yhi(1)
          ii = i*2
          
          capy(i,j) = half*(fapy(ii,jj)+fapy(ii+1,jj))
          if (capy(i,j) .ne. zero) then
             apinv = one/capy(i,j)
             cfcy(i,j) = half*apinv* &
                  ( fapy(ii  ,jj)*(half*ffcy(ii  ,jj)-fourth) &
                  + fapy(ii+1,jj)*(half*ffcy(ii+1,jj)+fourth) )
          end if
       end do
    end do

  end subroutine amrex_eb2_coarsen_from_fine


  subroutine amrex_eb2_build_cellflag_from_ap (lo, hi, &
       apx, xlo, xhi, apy, ylo, yhi, cflag, flo, fhi) &
       bind(c,name='amrex_eb2_build_cellflag_from_ap')
    integer, dimension(2), intent(in) :: lo, hi, xlo, xhi, ylo, yhi, flo, fhi
    real(amrex_real), intent(in   ) ::  apx (xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) ::  apy (ylo(1):yhi(1),ylo(2):yhi(2))
    integer         , intent(inout) :: cflag(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i,j, flg

    ! By default, all neighbors are already set.
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          flg = cflag(i,j)

          if (apx(i  ,j  ).eq.zero) flg = clear_neighbor(flg, -1,  0)
          if (apx(i+1,j  ).eq.zero) flg = clear_neighbor(flg,  1,  0)
          if (apy(i  ,j  ).eq.zero) flg = clear_neighbor(flg,  0, -1)
          if (apy(i  ,j+1).eq.zero) flg = clear_neighbor(flg,  0,  1)

          if (apx(i,j).ne.zero .and. apy(i-1,j).ne.zero) then
          else if (apx(i,j-1).ne.zero .and. apy(i,j).ne.zero) then
          else
             flg = clear_neighbor(flg,-1,-1)
          end if

          if (apx(i+1,j).ne.zero .and. apy(i+1,j).ne.zero) then
          else if (apx(i+1,j-1).ne.zero .and. apy(i,j).ne.zero) then
          else
             flg = clear_neighbor(flg,1,-1)
          end if

          if (apx(i,j).ne.zero .and. apy(i-1,j+1).ne.zero) then
          else if (apx(i,j+1).ne.zero .and. apy(i,j+1).ne.zero) then
          else
             flg = clear_neighbor(flg,-1,1)
          end if

          if (apx(i+1,j).ne.zero .and. apy(i+1,j+1).ne.zero) then
          else if (apx(i+1,j+1).ne.zero .and. apy(i,j+1).ne.zero) then
          else
             flg = clear_neighbor(flg,1,1)
          end if

          cflag(i,j) = flg
       end do
    end do
  end subroutine amrex_eb2_build_cellflag_from_ap

  
  subroutine amrex_eb2_check_mvmc (cclo, cchi, ndlo, ndhi, cls, clo, chi, fls, flo, fhi, ierr) &
       bind(c,name='amrex_eb2_check_mvmc')
    integer, dimension(2), intent(in) :: cclo, cchi, ndlo, ndhi, clo, chi, flo, fhi
    real(amrex_real), intent(inout) :: cls(clo(1):chi(1),clo(2):chi(2))
    real(amrex_real), intent(in   ) :: fls(flo(1):fhi(1),flo(2):fhi(2))
    integer, intent(inout) :: ierr

    integer :: i,j, ii, jj, ncuts

    do    j = ndlo(2), ndhi(2)
       do i = ndlo(1), ndhi(1)
          cls(i,j) = fls(2*i,2*j)
       end do
    end do

    do    j = cclo(2), cchi(2)
       jj = 2*j
       do i = cclo(1), cchi(1)
          ii = 2*i
          ncuts = 0
          if (has_cut(fls(ii  ,jj  ),fls(ii+1,jj  ))) ncuts = ncuts+1
          if (has_cut(fls(ii+1,jj  ),fls(ii+2,jj  ))) ncuts = ncuts+1
          if (has_cut(fls(ii  ,jj+2),fls(ii+1,jj+2))) ncuts = ncuts+1
          if (has_cut(fls(ii+1,jj+2),fls(ii+2,jj+2))) ncuts = ncuts+1
          if (has_cut(fls(ii  ,jj  ),fls(ii  ,jj+1))) ncuts = ncuts+1
          if (has_cut(fls(ii  ,jj+1),fls(ii  ,jj+2))) ncuts = ncuts+1
          if (has_cut(fls(ii+2,jj  ),fls(ii+2,jj+1))) ncuts = ncuts+1
          if (has_cut(fls(ii+2,jj+1),fls(ii+2,jj+2))) ncuts = ncuts+1
          if (ncuts .ne. 0 .and. ncuts .ne. 2) then
             ierr = 2
             return
          end if
       end do
    end do
  contains
    pure function has_cut(a,b)
      real(amrex_real), intent(in) :: a,b
      logical has_cut
      has_cut = (a.ge.zero .and. b.lt.zero) .or. (b.ge.zero .and. a.lt.zero)
    end function has_cut
  end subroutine amrex_eb2_check_mvmc

end module amrex_eb2_2d_moudle
