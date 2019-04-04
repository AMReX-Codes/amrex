
module amrex_eb2_3d_module

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module, only : zero, one, two, three, four, five, six, seven, eight,&
       nine, ten, eleven, fifteen, sixteen, half, third, fourth, sixth, eighth, twelfth
  use amrex_ebcellflag_module, only : regular_cell => regular, covered_cell => covered, &
       is_regular_cell, is_single_valued_cell, is_covered_cell, get_cell_type, get_neighbor_cells, &
       set_regular_cell, set_single_valued_cell, set_covered_cell, &
       set_neighbor, clear_neighbor, clear_allneighbors, is_neighbor
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
       fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, &
       ex, exlo, exhi, ey, eylo, eyhi, ez, ezlo, ezhi) &
       bind(c,name='amrex_eb2_gfab_build_types')
    integer, dimension(3), intent(in) :: lo, hi, slo, shi, clo, chi, &
         fxlo, fxhi, fylo, fyhi, fzlo, fzhi, exlo, exhi, eylo, eyhi, ezlo, ezhi
    real(amrex_real), intent(in   ) ::    s( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    integer(c_int)  , intent(inout) :: cell( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    integer(c_int)  , intent(inout) ::   fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
    integer(c_int)  , intent(inout) ::   fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
    integer(c_int)  , intent(inout) ::   fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
    integer(c_int)  , intent(inout) ::   ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3))
    integer(c_int)  , intent(inout) ::   ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3))
    integer(c_int)  , intent(inout) ::   ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3))

    integer :: i,j,k

    ! s > 0: body, s < 0: fluid

    do       k = lo(3)-2, hi(3)+2
       do    j = lo(2)-2, hi(2)+2
          do i = lo(1)-2, hi(1)+2
             if (       s(i,j  ,k  ).ge.zero .and. s(i+1,j  ,k  ).ge.zero &
                  .and. s(i,j+1,k  ).ge.zero .and. s(i+1,j+1,k  ).ge.zero &
                  .and. s(i,j  ,k+1).ge.zero .and. s(i+1,j  ,k+1).ge.zero &
                  .and. s(i,j+1,k+1).ge.zero .and. s(i+1,j+1,k+1).ge.zero ) then
                call set_covered_cell(cell(i,j,k))
             else if (  s(i,j  ,k  ).lt.zero .and. s(i+1,j  ,k  ).lt.zero &
                  .and. s(i,j+1,k  ).lt.zero .and. s(i+1,j+1,k  ).lt.zero &
                  .and. s(i,j  ,k+1).lt.zero .and. s(i+1,j  ,k+1).lt.zero &
                  .and. s(i,j+1,k+1).lt.zero .and. s(i+1,j+1,k+1).lt.zero ) then
                call set_regular_cell(cell(i,j,k))
             else
                call set_single_valued_cell(cell(i,j,k))
             end if
          end do
       end do
    end do

    ! x-face
    do       k = lo(3)-2, hi(3)+2
       do    j = lo(2)-2, hi(2)+2
          do i = lo(1)-2, hi(1)+3
             if (       s(i,j,k  ).ge.zero .and. s(i,j+1,k  ).ge.zero &
                  .and. s(i,j,k+1).ge.zero .and. s(i,j+1,k+1).ge.zero ) then
                fx(i,j,k) = covered
             else if (  s(i,j,k  ).lt.zero .and. s(i,j+1,k  ).lt.zero &
                  .and. s(i,j,k+1).lt.zero .and. s(i,j+1,k+1).lt.zero ) then
                fx(i,j,k) = regular
             else
                fx(i,j,k) = irregular
             end if
          end do
       end do
    end do

    ! y-face
    do       k = lo(3)-2, hi(3)+2
       do    j = lo(2)-2, hi(2)+3
          do i = lo(1)-2, hi(1)+2
             if (       s(i,j,k  ).ge.zero .and. s(i+1,j,k  ).ge.zero &
                  .and. s(i,j,k+1).ge.zero .and. s(i+1,j,k+1).ge.zero ) then
                fy(i,j,k) = covered
             else if (  s(i,j,k  ).lt.zero .and. s(i+1,j,k  ).lt.zero &
                  .and. s(i,j,k+1).lt.zero .and. s(i+1,j,k+1).lt.zero ) then
                fy(i,j,k) = regular
             else
                fy(i,j,k) = irregular
             end if
          end do
       end do
    end do

    ! z-face
    do       k = lo(3)-2, hi(3)+3
       do    j = lo(2)-2, hi(2)+2
          do i = lo(1)-2, hi(1)+2
             if (       s(i,j  ,k).ge.zero .and. s(i+1,j  ,k).ge.zero &
                  .and. s(i,j+1,k).ge.zero .and. s(i+1,j+1,k).ge.zero ) then
                fz(i,j,k) = covered
             else if (  s(i,j  ,k).lt.zero .and. s(i+1,j  ,k).lt.zero &
                  .and. s(i,j+1,k).lt.zero .and. s(i+1,j+1,k).lt.zero ) then
                fz(i,j,k) = regular
             else
                fz(i,j,k) = irregular
             end if
          end do
       end do
    end do

    ! x-edge
    do       k = lo(3)-2, hi(3)+3
       do    j = lo(2)-2, hi(2)+3
          do i = lo(1)-2, hi(1)+2
             if (s(i,j,k).ge.zero .and. s(i+1,j,k).ge.zero) then
                ex(i,j,k) = covered
             else if (s(i,j,k).lt.zero .and. s(i+1,j,k).lt.zero) then
                ex(i,j,k) = regular
             else
                ex(i,j,k) = irregular
             end if
          end do
       end do
    end do

    ! y-edge
    do       k = lo(3)-2, hi(3)+3
       do    j = lo(2)-2, hi(2)+2
          do i = lo(1)-2, hi(1)+3
             if (s(i,j,k).ge.zero .and. s(i,j+1,k).ge.zero) then
                ey(i,j,k) = covered
             else if (s(i,j,k).lt.zero .and. s(i,j+1,k).lt.zero) then
                ey(i,j,k) = regular
             else
                ey(i,j,k) = irregular
             end if
          end do
       end do
    end do

    ! z-edge
    do       k = lo(3)-2, hi(3)+2
       do    j = lo(2)-2, hi(2)+3
          do i = lo(1)-2, hi(1)+3
             if (s(i,j,k).ge.zero .and. s(i,j,k+1).ge.zero) then
                ez(i,j,k) = covered
             else if (s(i,j,k).lt.zero .and. s(i,j,k+1).lt.zero) then
                ez(i,j,k) = regular
             else
                ez(i,j,k) = irregular
             end if
          end do
       end do
    end do

  end subroutine amrex_eb2_gfab_build_types


  subroutine amrex_eb2_build_faces (lo, hi, cell, clo, chi, &
       fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, &
       ex, exlo, exhi, ey, eylo, eyhi, ez, ezlo, ezhi, &
       levset, slo, shi,&
       interx, ixlo, ixhi, intery, iylo, iyhi, interz, izlo, izhi, &
       apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi, &
       fcx, cxlo, cxhi, fcy, cylo, cyhi, fcz, czlo, czhi, &
       m2x, mxlo, mxhi, m2y, mylo, myhi, m2z, mzlo, mzhi, &
       dx, dxinv, problo) &
       bind(c,name='amrex_eb2_build_faces')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, slo, shi, &
         fxlo, fxhi, fylo, fyhi, fzlo, fzhi, exlo, exhi, eylo, eyhi, ezlo, ezhi, &
         ixlo, ixhi, iylo, iyhi, izlo, izhi, axlo, axhi, aylo, ayhi, azlo, azhi, &
         cxlo, cxhi, cylo, cyhi, czlo, czhi, mxlo, mxhi, mylo, myhi, mzlo, mzhi
    real(amrex_real), dimension(3) :: dx, dxinv, problo
    integer(c_int)  , intent(inout) ::   cell( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    integer(c_int)  , intent(inout) ::     fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
    integer(c_int)  , intent(inout) ::     fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
    integer(c_int)  , intent(inout) ::     fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
    integer(c_int)  , intent(in   ) ::     ex(exlo(1):exhi(1),exlo(2):exhi(2),exlo(3):exhi(3))
    integer(c_int)  , intent(in   ) ::     ey(eylo(1):eyhi(1),eylo(2):eyhi(2),eylo(3):eyhi(3))
    integer(c_int)  , intent(in   ) ::     ez(ezlo(1):ezhi(1),ezlo(2):ezhi(2),ezlo(3):ezhi(3))
    real(amrex_real), intent(in   ) :: levset( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    real(amrex_real), intent(in   ) :: interx(ixlo(1):ixhi(1),ixlo(2):ixhi(2),ixlo(3):ixhi(3))
    real(amrex_real), intent(in   ) :: intery(iylo(1):iyhi(1),iylo(2):iyhi(2),iylo(3):iyhi(3))
    real(amrex_real), intent(in   ) :: interz(izlo(1):izhi(1),izlo(2):izhi(2),izlo(3):izhi(3))
    real(amrex_real), intent(inout) ::    apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(amrex_real), intent(inout) ::    apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(amrex_real), intent(inout) ::    apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(amrex_real), intent(inout) ::    fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
    real(amrex_real), intent(inout) ::    fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2)
    real(amrex_real), intent(inout) ::    fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2)
    real(amrex_real), intent(inout) ::    m2x(mxlo(1):mxhi(1),mxlo(2):mxhi(2),mxlo(3):mxhi(3),3)
    real(amrex_real), intent(inout) ::    m2y(mylo(1):myhi(1),mylo(2):myhi(2),mylo(3):myhi(3),3)
    real(amrex_real), intent(inout) ::    m2z(mzlo(1):mzhi(1),mzlo(2):mzhi(2),mzlo(3):mzhi(3),3)

    integer :: i,j,k, ncuts
    real(amrex_real) :: cut
    real(amrex_real) :: lxm, lxp, lym, lyp, lzm, lzp, bcx, bcy, bcz
    real(amrex_real), parameter :: sixteenth = 1.d0/16.d0
    real(amrex_real), parameter :: twentyfourth = 1.d0/24.d0

    ! x-face
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+2
             if (fx(i,j,k) .eq. regular) then
!                apx(i,j,k) = one
!                fcx(i,j,k,:) = zero
                m2x(i,j,k,1) = twelfth
                m2x(i,j,k,2) = twelfth
                m2x(i,j,k,3) = zero
             else if (fx(i,j,k) .eq. covered) then
                apx(i,j,k) = zero
                fcx(i,j,k,:) = zero
                m2x(i,j,k,:) = zero
             else
                ncuts = 0
                bcy = zero
                bcz = zero

                if (ey(i,j,k) .eq. regular) then
                   lym = one
                else if (ey(i,j,k) .eq. covered) then
                   lym = zero
                else
                   ncuts = ncuts+1
                   cut = (intery(i,j,k)-(problo(2)+j*dx(2)))*dxinv(2)
                   bcy = bcy + cut
                   if (levset(i,j,k) .lt. zero) then
                      lym = cut
                   else
                      lym = one-cut
                   end if
                   lym = min(max(zero,lym),one)
                end if

                if (ey(i,j,k+1) .eq. regular) then
                   lyp = one
                else if (ey(i,j,k+1) .eq. covered) then
                   lyp = zero
                else
                   ncuts = ncuts+1
                   cut = (intery(i,j,k+1)-(problo(2)+j*dx(2)))*dxinv(2)
                   bcy = bcy + cut
                   bcz = bcz + one
                   if (levset(i,j,k+1) .lt. zero) then
                      lyp = cut
                   else
                      lyp = one-cut
                   end if
                   lyp = min(max(zero,lyp),one)
                end if

                if (ez(i,j,k) .eq. regular) then
                   lzm = one
                else if (ez(i,j,k) .eq. covered) then
                   lzm = zero
                else
                   ncuts = ncuts+1
                   cut = (interz(i,j,k)-(problo(3)+k*dx(3)))*dxinv(3)
                   bcz = bcz + cut
                   if (levset(i,j,k) .lt. zero) then
                      lzm = cut
                   else
                      lzm = one-cut
                   end if
                   lzm = min(max(zero,lzm),one)
                end if

                if (ez(i,j+1,k) .eq. regular) then
                   lzp = one
                else if (ez(i,j+1,k) .eq. covered) then
                   lzp = zero
                else
                   ncuts = ncuts+1
                   cut = (interz(i,j+1,k)-(problo(3)+k*dx(3)))*dxinv(3)
                   bcy = bcy + one
                   bcz = bcz + cut
                   if (levset(i,j+1,k) .lt. zero) then
                      lzp = cut
                   else
                      lzp = one-cut
                   end if
                end if

                if (ncuts .gt. 2) then
#ifdef AMREX_DEBUG
                   print *, "ncuts = ", i,j,k,ncuts
                   flush(6)
#endif
                   call amrex_error("amrex_eb2_build_faces: more than 2 cuts not suported")
                end if

                if (lym.le.small .and. lyp.le.small .and. lzm.le.small .and. lzp.le.small) then
                   apx(i,j,k) = zero
                   fcx(i,j,k,:) = zero
                   m2x(i,j,k,1) = zero
                   m2x(i,j,k,2) = zero
                   m2x(i,j,k,3) = zero
                else if (lym.eq.lyp .and. lzm.eq.lzp) then
                   apx(i,j,k) = one
                   fcx(i,j,k,:) = zero
                   m2x(i,j,k,1) = twelfth
                   m2x(i,j,k,2) = twelfth
                   m2x(i,j,k,3) = zero
                else
                   bcy = half*bcy - half
                   bcz = half*bcz - half
                   call cut_face_2d(apx(i,j,k),fcx(i,j,k,1),fcx(i,j,k,2), &
                        m2x(i,j,k,1),m2x(i,j,k,2),m2x(i,j,k,3),lzm,lzp,lym,lyp,bcy,bcz)
                end if

                if (apx(i,j,k) .eq. zero) then
                   fx(i,j,k) = covered
                else if (apx(i,j,k) .eq. one) then
                   fx(i,j,k) = regular
                end if

             end if
          end do
       end do
    end do

    ! y-face
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+2
          do i = lo(1)-1, hi(1)+1
             if (fy(i,j,k) .eq. regular) then
!                apy(i,j,k) = one
!                fcy(i,j,k,:) = zero
                m2y(i,j,k,1) = twelfth
                m2y(i,j,k,2) = twelfth
                m2y(i,j,k,3) = zero
             else if (fy(i,j,k) .eq. covered) then
                apy(i,j,k) = zero
                fcy(i,j,k,:) = zero
                m2y(i,j,k,:) = zero
             else
                ncuts = 0
                bcx = zero
                bcz = zero

                if (ex(i,j,k) .eq. regular) then
                   lxm = one
                else if (ex(i,j,k) .eq. covered) then
                   lxm = zero
                else
                   ncuts = ncuts+1
                   cut = (interx(i,j,k)-(problo(1)+i*dx(1)))*dxinv(1)
                   bcx = bcx + cut
                   if (levset(i,j,k) .lt. zero) then
                      lxm = cut
                   else
                      lxm = one-cut
                   end if
                   lxm = min(max(zero,lxm),one)
                end if

                if (ex(i,j,k+1) .eq. regular) then
                   lxp = one
                else if (ex(i,j,k+1) .eq. covered) then
                   lxp = zero
                else
                   ncuts = ncuts+1
                   cut = (interx(i,j,k+1)-(problo(1)+i*dx(1)))*dxinv(1)
                   bcx = bcx + cut
                   bcz = bcz + one
                   if (levset(i,j,k+1) .lt. zero) then
                      lxp = cut
                   else
                      lxp = one-cut
                   end if
                   lxp = min(max(zero,lxp),one)
                end if

                if (ez(i,j,k) .eq. regular) then
                   lzm = one
                else if (ez(i,j,k) .eq. covered) then
                   lzm = zero
                else
                   ncuts = ncuts+1
                   cut = (interz(i,j,k)-(problo(3)+k*dx(3)))*dxinv(3)
                   bcz = bcz + cut
                   if (levset(i,j,k) .lt. zero) then
                      lzm = cut
                   else
                      lzm = one-cut
                   end if
                end if

                if (ez(i+1,j,k) .eq. regular) then
                   lzp = one
                else if (ez(i+1,j,k) .eq. covered) then
                   lzp = zero
                else
                   ncuts = ncuts+1
                   cut = (interz(i+1,j,k)-(problo(3)+k*dx(3)))*dxinv(3)
                   bcx = bcx + one
                   bcz = bcz + cut
                   if (levset(i+1,j,k) .lt. zero) then
                      lzp = cut
                   else
                      lzp = one-cut
                   end if
                end if

                if (ncuts .gt. 2) then
#ifdef AMREX_DEBUG
                   print *, "ncuts = ", i,j,k,ncuts
                   flush(6)
#endif
                   call amrex_error("amrex_eb2_build_faces: more than 2 cuts not supported")
                end if

                if (lxm.le.small .and. lxp.le.small .and. lzm.le.small .and. lzp.le.small) then
                   apy(i,j,k) = zero
                   fcy(i,j,k,:) = zero
                   m2y(i,j,k,1) = zero
                   m2y(i,j,k,2) = zero
                   m2y(i,j,k,3) = zero
                else if (lxm.eq.lxp .and. lzm.eq.lzp) then
                   apy(i,j,k) = one
                   fcy(i,j,k,:) = zero
                   m2y(i,j,k,1) = twelfth
                   m2y(i,j,k,2) = twelfth
                   m2y(i,j,k,3) = zero
                else
                   bcx = half*bcx - half
                   bcz = half*bcz - half
                   call cut_face_2d(apy(i,j,k),fcy(i,j,k,1),fcy(i,j,k,2),&
                        m2y(i,j,k,1),m2y(i,j,k,2),m2y(i,j,k,3),lzm,lzp,lxm,lxp,bcx,bcz)
                end if

                if (apy(i,j,k) .eq. zero) then
                   fy(i,j,k) = covered
                else if (apy(i,j,k) .eq. one) then
                   fy(i,j,k) = regular
                end if

             end if
          end do
       end do
    end do

    ! z-face
    do       k = lo(3)-1, hi(3)+2
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             if (fz(i,j,k) .eq. regular) then
!                apz(i,j,k) = one
!                fcz(i,j,k,:) = zero
                m2z(i,j,k,1) = twelfth
                m2z(i,j,k,2) = twelfth
                m2z(i,j,k,3) = zero
             else if (fz(i,j,k) .eq. covered) then
                apz(i,j,k) = zero
                fcz(i,j,k,:) = zero
                m2z(i,j,k,:) = zero
             else

                ncuts = 0
                bcx = zero
                bcy = zero

                if (ex(i,j,k) .eq. regular) then
                   lxm = one
                else if (ex(i,j,k) .eq. covered) then
                   lxm = zero
                else
                   ncuts = ncuts+1
                   cut = (interx(i,j,k)-(problo(1)+i*dx(1)))*dxinv(1)
                   bcx = bcx + cut
                   if (levset(i,j,k) .lt. zero) then
                      lxm = cut
                   else
                      lxm = one-cut
                   end if
                   lxm = min(max(zero,lxm),one)
                end if

                if (ex(i,j+1,k) .eq. regular) then
                   lxp = one
                else if (ex(i,j+1,k) .eq. covered) then
                   lxp = zero
                else
                   ncuts = ncuts+1
                   cut = (interx(i,j+1,k)-(problo(1)+i*dx(1)))*dxinv(1)
                   bcx = bcx + cut
                   bcy = bcy + one
                   if (levset(i,j+1,k) .lt. zero) then
                      lxp = cut
                   else
                      lxp = one-cut
                   end if
                   lxp = min(max(zero,lxp),one)
                end if

                if (ey(i,j,k) .eq. regular) then
                   lym = one
                else if (ey(i,j,k) .eq. covered) then
                   lym = zero
                else
                   ncuts = ncuts+1
                   cut = (intery(i,j,k)-(problo(2)+j*dx(2)))*dxinv(2)
                   bcy = bcy + cut
                   if (levset(i,j,k) .lt. zero) then
                      lym = cut
                   else
                      lym = one-cut
                   end if
                end if

                if (ey(i+1,j,k) .eq. regular) then
                   lyp = one
                else if (ey(i+1,j,k) .eq. covered) then
                   lyp = zero
                else
                   ncuts = ncuts+1
                   cut = (intery(i+1,j,k)-(problo(2)+j*dx(2)))*dxinv(2)
                   bcx = bcx + one
                   bcy = bcy + cut
                   if (levset(i+1,j,k) .lt. zero) then
                      lyp = cut
                   else
                      lyp = one-cut
                   end if
                end if

                if (ncuts .gt. 2) then
#ifdef AMREX_DEBUG
                   print *, "ncuts = ", i,j,k,ncuts
                   flush(6)
#endif
                   call amrex_error("amrex_eb2_build_faces: more than 2 cuts not supported")
                end if

                if (lxm.le.small .and. lxp.le.small .and. lym.le.small .and. lyp.le.small) then
                   apz(i,j,k) = zero
                   fcz(i,j,k,:) = zero
                   m2z(i,j,k,1) = zero
                   m2z(i,j,k,2) = zero
                   m2z(i,j,k,3) = zero
                else if (lxm.eq.lxp .and. lym.eq.lyp) then
                   apz(i,j,k) = one
                   fcz(i,j,k,:) = zero
                   m2z(i,j,k,1) = twelfth
                   m2z(i,j,k,2) = twelfth
                   m2z(i,j,k,3) = zero
                else
                   bcx = half*bcx - half
                   bcy = half*bcy - half
                   call cut_face_2d(apz(i,j,k),fcz(i,j,k,1),fcz(i,j,k,2), &
                        m2z(i,j,k,1),m2z(i,j,k,2),m2z(i,j,k,3),lym,lyp,lxm,lxp,bcx,bcy)
                end if

                if (apz(i,j,k) .eq. zero) then
                   fz(i,j,k) = covered
                else if (apz(i,j,k) .eq. one) then
                   fz(i,j,k) = regular
                end if

             end if
          end do
       end do
    end do

    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             if (is_single_valued_cell(cell(i,j,k))) then
                if  (fx(i,j,k) .eq. covered .and. fx(i+1,j,k) .eq. covered .and. &
                     fy(i,j,k) .eq. covered .and. fy(i,j+1,k) .eq. covered .and. &
                     fz(i,j,k) .eq. covered .and. fz(i,j,k+1) .eq. covered ) then
                   call set_covered_cell(cell(i,j,k))
                else if (fx(i,j,k) .eq. regular .and. fx(i+1,j,k) .eq. regular .and. &
                     &   fy(i,j,k) .eq. regular .and. fy(i,j+1,k) .eq. regular .and. &
                     &   fz(i,j,k) .eq. regular .and. fz(i,j,k+1) .eq. regular ) then
                   call set_regular_cell(cell(i,j,k))
                end if
             end if
          end do
       end do
    end do

  contains

    subroutine cut_face_2d (areafrac, centx, centy, Sx2, Sy2, Sxy, axm, axp, aym, ayp, bcx, bcy)
      real(amrex_real), intent(out) :: areafrac, centx, centy, Sx2, Sy2, Sxy
      real(amrex_real), intent(in) :: axm, axp, aym, ayp, bcx, bcy

      real(amrex_real) :: apnorm, apnorminv, nx, ny, x_ym, x_yp, y_xm, y_xp, aa, af1, af2
      real(amrex_real) :: dx, dx2, dx3, dx4, dy, dy2, dy3, dy4, nxabs, nyabs, S_b, signx, signy
      real(amrex_real), parameter :: tiny = 1.d-15
      real(amrex_real), parameter :: almostone = 1.d0-1.d-15

      apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2)
      apnorminv = one/apnorm
      nx = (axm-axp) * apnorminv  ! pointing to the wall
      ny = (aym-ayp) * apnorminv

      nxabs = abs(nx)
      nyabs = abs(ny)

      if (nxabs .lt. tiny .or. nyabs .gt. almostone) then
         areafrac = half*(axm+axp)
         if (areafrac .gt. one-small) then
            areafrac = one
            centx = zero
            centy = zero
            Sx2 = twelfth
            Sy2 = twelfth
            Sxy = zero
         else if (areafrac .lt. small) then
            areafrac = zero
            centx = zero
            centy = zero
            Sx2 = zero
            Sy2 = zero
            Sxy = zero
         else
            centx = zero
            centy = (eighth*(ayp-aym) + ny*half*bcy**2)/areafrac
            Sx2 = twentyfourth*(axm+axp)
            Sy2 = twentyfourth*(ayp+aym) + ny*third*bcy**3
            Sxy = zero
         end if
      else if (nyabs .lt. tiny .or. nxabs .gt. almostone) then
         areafrac = half*(aym+ayp)
         if (areafrac .gt. one-small) then
            areafrac = one
            centx = zero
            centy = zero
            Sx2 = twelfth
            Sy2 = twelfth
            Sxy = zero
         else if (areafrac .lt. small) then
            areafrac = zero
            centx = zero
            centy = zero
            Sx2 = zero
            Sy2 = zero
            Sxy = zero
         else
            centx = (eighth*(axp-axm) + nx*half*bcx**2)/areafrac
            centy = zero
            Sx2 = twentyfourth*(axp+axm) + nx*third*bcx**3
            Sy2 = twentyfourth*(ayp+aym)
            Sxy = zero
         end if
      else
         if (nx .gt. zero) then
            x_ym = -half + aym
            x_yp = -half + ayp
            signx = one
         else
            x_ym = half - aym
            x_yp = half - ayp
            signx = -one
         end if
         aa = nxabs/ny
         dx  = x_ym - x_yp
         dx2 = dx * (x_ym + x_yp)
         dx3 = dx * (x_ym**2 + x_ym*x_yp + x_yp**2)
         dx4 = dx * (x_ym + x_yp) * (x_ym**2 + x_yp**2)
         af1 = half*(axm+axp) + aa*half*dx2
         centx = eighth*(axp-axm) + aa*sixth*dx3
         Sx2 = twentyfourth*(axm+axp) + aa*twelfth*dx4

         if (ny .gt. zero) then
            y_xm = -half + axm
            y_xp = -half + axp
            signy = one
         else
            y_xm = half - axm
            y_xp = half - axp
            signy = -one
         end if
         aa = nyabs/nx
         dy = y_xm - y_xp
         dy2 = dy * (y_xm + y_xp)
         dy3 = dy * (y_xm**2 + y_xm*y_xp + y_xp**2)
         dy4 = dy * (y_xm + y_xp) * (y_xm**2 + y_xp**2)
         af2 = half*(aym+ayp) + aa*half*dy2
         centy = eighth*(ayp-aym) + aa*sixth*dy3
         Sy2 = twentyfourth*(aym+ayp) + aa*twelfth*dy4

         if (nxabs .lt. nyabs) then
            S_b = (Sx2 - twentyfourth - signx*sixth*(x_ym**3+x_yp**3)) / ny
            Sxy = -signy*sixteenth*dy2 + half*nx*S_b
         else
            S_b = (Sy2 - twentyfourth - signy*sixth*(y_xm**3+y_xp**3)) / nx
            Sxy = -signx*sixteenth*dx2 + half*ny*S_b
         end if

         areafrac = half*(af1+af2)
         if (areafrac .gt. one-small) then
            areafrac = one
            centx = zero
            centy = zero
            Sx2 = twelfth
            Sy2 = twelfth
            Sxy = zero
         else if (areafrac .lt. small) then
            areafrac = zero
            centx = zero
            centy = zero
            Sx2 = zero
            Sy2 = zero
            Sxy = zero
         else
            centx = centx*(one/areafrac)
            centy = centy*(one/areafrac)
            centx = min(max(centx,-half),half)
            centy = min(max(centy,-half),half)
         end if
      end if

    end subroutine cut_face_2d

  end subroutine amrex_eb2_build_faces


  subroutine amrex_eb2_build_cells (lo, hi, cell, clo, chi, &
       fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, &
       apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi, &
       fcx, cxlo, cxhi, fcy, cylo, cyhi, fcz, czlo, czhi, &
       m2x, mxlo, mxhi, m2y, mylo, myhi, m2z, mzlo, mzhi, &
       vfrac, vlo, vhi, vcent, tlo, thi, barea, alo, ahi, &
       bcent, blo, bhi, bnorm, mlo, mhi) &
       bind(c,name='amrex_eb2_build_cells')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, &
         fxlo, fxhi, fylo, fyhi, fzlo, fzhi, axlo, axhi, aylo, ayhi, azlo, azhi, &
         cxlo, cxhi, cylo, cyhi, czlo, czhi, mxlo, mxhi, mylo, myhi, mzlo, mzhi, &
         vlo, vhi, tlo, thi, alo, ahi, blo, bhi, mlo, mhi
    integer(c_int)  , intent(inout) ::  cell( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    integer(c_int)  , intent(inout) ::    fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
    integer(c_int)  , intent(inout) ::    fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
    integer(c_int)  , intent(inout) ::    fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
    real(amrex_real), intent(inout) ::   apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(amrex_real), intent(inout) ::   apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(amrex_real), intent(inout) ::   apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(amrex_real), intent(in   ) ::   fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
    real(amrex_real), intent(in   ) ::   fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2)
    real(amrex_real), intent(in   ) ::   fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2)
    real(amrex_real), intent(in   ) ::   m2x(mxlo(1):mxhi(1),mxlo(2):mxhi(2),mxlo(3):mxhi(3),3)
    real(amrex_real), intent(in   ) ::   m2y(mylo(1):myhi(1),mylo(2):myhi(2),mylo(3):myhi(3),3)
    real(amrex_real), intent(in   ) ::   m2z(mzlo(1):mzhi(1),mzlo(2):mzhi(2),mzlo(3):mzhi(3),3)
    real(amrex_real), intent(inout) :: vfrac( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3))
    real(amrex_real), intent(inout) :: vcent( tlo(1): thi(1), tlo(2): thi(2), tlo(3): thi(3),3)
    real(amrex_real), intent(inout) :: barea( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
    real(amrex_real), intent(inout) :: bcent( blo(1): bhi(1), blo(2): bhi(2), blo(3): bhi(3),3)
    real(amrex_real), intent(inout) :: bnorm( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3),3)

    integer :: i,j,k, ngbr(-1:1,-1:1,-1:1), flg
    real(amrex_real) :: axm, axp, aym, ayp, azm, azp, aax, aay, aaz, B0, Bx, By, Bz
    real(amrex_real) :: dapx, dapy, dapz, apnorm, apnorminv, nx, ny, nz, bainv
    real(amrex_real) :: ny2, ny3, ny4, nz2, nz3, nz4, nz5, Sx, Sy, Sz, den
    real(amrex_real) :: rhs(18), b(9)

    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             if (is_regular_cell(cell(i,j,k))) then
!                vfrac(i,j,k) = one
!                vcent(i,j,k,:) = zero
!                bcent(i,j,k,:) = -one
!                bnorm(i,j,k,:) = zero
!                barea(i,j,k) = zero
             else if (is_covered_cell(cell(i,j,k))) then
                vfrac(i,j,k) = zero
!                vcent(i,j,k,:) = zero
!                bcent(i,j,k,:) = -one
!                bnorm(i,j,k,:) = zero
!                barea(i,j,k) = zero
             else

                axm = apx(i,j,k)
                axp = apx(i+1,j,k)
                aym = apy(i,j,k)
                ayp = apy(i,j+1,k)
                azm = apz(i,j,k)
                azp = apz(i,j,k+1)

                ! Check for multple cuts
                ! We know there are no multiple cuts on faces by now.
                ! So we only need to check the case that there are two cuts
                ! at the opposite corners.
                if ( axm .ge. half .and. axm .lt. one .and. &
                     axp .ge. half .and. axp .lt. one .and. &
                     aym .ge. half .and. aym .lt. one .and. &
                     ayp .ge. half .and. ayp .lt. one .and. &
                     azm .ge. half .and. azm .lt. one .and. &
                     azp .ge. half .and. azp .lt. one ) then
#ifdef AMREX_DEBUG
                   print *, "amrex_eb2_build_cells: multiple cuts in cell ", i,j,k
#endif
                   call amrex_error("amrex_eb2_build_cells: multiple cuts not supported")
                end if

                dapx = axm - axp
                dapy = aym - ayp
                dapz = azm - azp
                apnorm = sqrt(dapx**2 + dapy**2 + dapz**2)
                if (apnorm .eq. zero) then
#ifdef AMREX_DEBUG
                   print *, "amrex_eb2_build_cells: multiple cuts in cell ", i,j,k
#endif
                   call amrex_error("amrex_eb2_build_cells: apnorm==0, multiple cuts not supported")
                end if
                apnorminv = one/apnorm
                nx = dapx * apnorminv
                ny = dapy * apnorminv
                nz = dapz * apnorminv
                bnorm(i,j,k,1) = nx
                bnorm(i,j,k,2) = ny
                bnorm(i,j,k,3) = nz
                barea(i,j,k) = nx*dapx + ny*dapy + nz*dapz

                aax = half*(axm+axp)
                aay = half*(aym+ayp)
                aaz = half*(azm+azp)
                B0 = aax + aay + aaz
                Bx = -nx*aax + ny*(aym*fcy(i,j,k,1)-ayp*fcy(i,j+1,k,1)) &
                     &       + nz*(azm*fcz(i,j,k,1)-azp*fcz(i,j,k+1,1))
                By = -ny*aay + nx*(axm*fcx(i,j,k,1)-axp*fcx(i+1,j,k,1)) &
                     &       + nz*(azm*fcz(i,j,k,2)-azp*fcz(i,j,k+1,2))
                Bz = -nz*aaz + nx*(axm*fcx(i,j,k,2)-axp*fcx(i+1,j,k,2)) &
                     &       + ny*(aym*fcy(i,j,k,2)-ayp*fcy(i,j+1,k,2))

                vfrac(i,j,k) = half*(B0 + nx*Bx + ny*By + nz*Bz)

                bainv = one / barea(i,j,k)
                bcent(i,j,k,1) = bainv * (Bx + nx*vfrac(i,j,k))
                bcent(i,j,k,2) = bainv * (By + ny*vfrac(i,j,k))
                bcent(i,j,k,3) = bainv * (Bz + nz*vfrac(i,j,k))

                rhs(1) = fourth*(axp-axm)
                rhs(2) = m2x(i+1,j,k,1) - m2x(i,j,k,1)
                rhs(3) = m2x(i+1,j,k,2) - m2x(i,j,k,2)
                rhs(4) = half*(axp*fcx(i+1,j,k,1) + axm*fcx(i,j,k,1))
                rhs(5) = half*(axp*fcx(i+1,j,k,2) + axm*fcx(i,j,k,2))
                rhs(6) = m2x(i+1,j,k,3) - m2x(i,j,k,3)
                !
                rhs(7) = m2y(i,j+1,k,1) - m2y(i,j,k,1)
                rhs(8) = fourth*(ayp-aym)
                rhs(9) = m2y(i,j+1,k,2) - m2y(i,j,k,2)
                rhs(10) = half*(ayp*fcy(i,j+1,k,1) + aym*fcy(i,j,k,1))
                rhs(11) = m2y(i,j+1,k,3) - m2y(i,j,k,3)
                rhs(12) = half*(ayp*fcy(i,j+1,k,2) + aym*fcy(i,j,k,2))
                !
                rhs(13) = m2z(i,j,k+1,1) - m2z(i,j,k,1)
                rhs(14) = m2z(i,j,k+1,2) - m2z(i,j,k,2)
                rhs(15) = fourth*(azp-azm)
                rhs(16) = m2z(i,j,k+1,3) - m2z(i,j,k,3)
                rhs(17) = half*(azp*fcz(i,j,k+1,1) + azm*fcz(i,j,k,1))
                rhs(18) = half*(azp*fcz(i,j,k+1,2) + azm*fcz(i,j,k,2))

                b(1) = two*rhs(1) + rhs(10) + rhs(17)
                b(2) = rhs(4) + two*rhs(8) + rhs(18)
                b(3) = rhs(5) + rhs(12) + two*rhs(15)
                b(4) = -nx*rhs(1) - ny*rhs(7) - nz*rhs(13)
                b(5) = -nx*rhs(2) - ny*rhs(8) - nz*rhs(14)
                b(6) = -nx*rhs(3) - ny*rhs(9) - nz*rhs(15)
                b(7) = -nx*rhs(4) - ny*rhs(10) - nz*rhs(16)
                b(8) = -nx*rhs(5) - ny*rhs(11) - nz*rhs(17)
                b(9) = -nx*rhs(6) - ny*rhs(12) - nz*rhs(18)

                ny2 = ny*ny
                ny3 = ny2*ny
                ny4 = ny3*ny
                nz2 = nz*nz
                nz3 = nz2*nz
                nz4 = nz3*nz
                nz5 = nz4*nz

                Sx = (five*(b(1)*(five - three*ny2) + two*b(4)*nx*(five - three*ny2) + &
                     ny*(nx*(b(2) + two*b(5)*ny) + b(7)*(six - four*ny2))) + &
                     (two*b(8)*(fifteen - eleven*ny2 + ny4) + &
                     nx*(b(3)*(five - two*ny2) - two*b(9)*ny*(-five + ny2)))*nz + &
                     (-22.d0*b(7)*ny - two*nx*(fifteen*b(4) - five*b(6) + b(2)*ny) + &
                     ny2*((sixteen*b(4) - four*(b(5) + b(6)))*nx + ten*b(7)*ny) + &
                     b(1)*(-fifteen + eight*ny2))*nz2 + &
                     two*(-(b(9)*nx*ny) + five*b(8)*(-two + ny2))*nz3 + &
                     two*b(7)*ny*nz4)

                Sy = (five*(two*b(7)*nx*(one + two*ny2) + b(2)*(two + three*ny2) + &
                     ny*(b(1)*nx - two*b(4)*(-one + ny2) + b(5)*(four + six*ny2))) + &
                     (two*b(9)*(five + nine*ny2 + ny4) + &
                     ny*(two*b(8)*nx*(four + ny2) + b(3)*(three + two*ny2)))*nz + &
                     (two*b(7)*nx*(four - five*ny2) - eight*b(2)*(-one + ny2) + &
                     two*ny*(-seven*b(4) + eight*b(5) + three*b(6) - b(1)*nx + &
                     two*(b(4) - four*b(5) + b(6))*ny2))*nz2 + &
                     two*(b(3)*ny + b(9)*(four - three*ny2))*nz3 + &
                     (-eight*(b(2) + b(7)*nx) + four*(b(4) - four*b(5) + b(6))*ny)*nz4 - &
                     eight*b(9)*nz5)

                Sz = (-two*(b(3) + b(8)*nx + b(9)*ny)*(-five - four*ny2 + four*ny4) + &
                     (five*(two*b(4) + four*b(6) + b(1)*nx) + (three*b(2) + eight*b(7)*nx)*ny - &
                     two*(seven*b(4) - three*b(5) - eight*b(6) + b(1)*nx)*ny2 + &
                     two*b(2)*ny3 + four*(b(4) + b(5) - four*b(6))*ny4)*nz + &
                     (b(3)*(fifteen - eight*ny2) - six*b(9)*ny*(-three + ny2) - &
                     ten*b(8)*nx*(-two + ny2))*nz2 + &
                     two*(-five*b(4) + fifteen*b(6) + (b(2) + b(7)*nx)*ny + &
                     two*(b(4) + b(5) - four*b(6))*ny2)*nz3 + two*b(9)*ny*nz4)

                den = one / (ten*(five + four*nz2 - four*nz4 + two*ny4*(-two + nz2) + &
                     two*ny2*(two - three*nz2 + nz4)) * (vfrac(i,j,k)+1.d-50) )

                vcent(i,j,k,1) = sx * den
                vcent(i,j,k,2) = Sy * den
                vcent(i,j,k,3) = Sz * den

                ! remove small cells
                if (vfrac(i,j,k) < small) then
                   vfrac(i,j,k) = zero
                   vcent(i,j,k,:) = zero
                   bcent(i,j,k,:) = -one
                   bnorm(i,j,k,:) = zero
                   barea(i,j,k) = zero
                   call set_covered_cell(cell(i,j,k))
                end if
             end if
          end do
       end do
    end do

    ! fix faces for small cells
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)  , hi(1)+1
             if (vfrac(i-1,j,k) < small .or. vfrac(i,j,k) < small) then
                fx(i,j,k) = covered
                apx(i,j,k) = zero
             end if
          end do
       end do
    end do
    !
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)  , hi(2)+1
          do i = lo(1)-1, hi(1)+1
             if (vfrac(i,j-1,k) < small .or. vfrac(i,j,k) < small) then
                fy(i,j,k) = covered
                apy(i,j,k) = zero
             end if
          end do
       end do
    end do
    !
    do       k = lo(3)  , hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             if (vfrac(i,j,k-1) < small .or. vfrac(i,j,k) < small) then
                fz(i,j,k) = covered
                apz(i,j,k) = zero
             end if
          end do
       end do
    end do

    ! Build neighbors.  By default, all 26 neighbors are already set.

    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             if (is_covered_cell(cell(i,j,k))) then
                cell(i,j,k) = clear_allneighbors(cell(i,j,k))
             else
                flg = cell(i,j,k)

                if (fx(i,j,k).eq.covered) then
                   flg = clear_neighbor(flg,-1,0,0)
                end if
                if (fx(i+1,j,k).eq.covered) then
                   flg = clear_neighbor(flg,1,0,0)
                end if
                if (fy(i,j,k).eq.covered) then
                   flg = clear_neighbor(flg,0,-1,0)
                end if
                if (fy(i,j+1,k).eq.covered) then
                   flg = clear_neighbor(flg,0,1,0)
                end if
                if (fz(i,j,k).eq.covered) then
                   flg = clear_neighbor(flg,0,0,-1)
                end if
                if (fz(i,j,k+1).eq.covered) then
                   flg = clear_neighbor(flg,0,0,1)
                end if
                
                ! x-y
                if (fx(i,j,k).ne.covered .and. fy(i-1,j,k).ne.covered) then
                else if (fx(i,j-1,k).ne.covered .and. fy(i,j,k).ne.covered) then
                else
                   flg = clear_neighbor(flg,-1,-1,0)
                end if
                
                if (fx(i+1,j,k).ne.covered .and. fy(i+1,j,k).ne.covered) then
                else if (fx(i+1,j-1,k).ne.covered .and. fy(i,j,k).ne.covered) then
                else
                   flg = clear_neighbor(flg,1,-1,0)
                end if
                
                if (fx(i,j,k).ne.covered .and. fy(i-1,j+1,k).ne.covered) then
                else if (fx(i,j+1,k).ne.covered .and. fy(i,j+1,k).ne.covered) then
                else
                   flg = clear_neighbor(flg,-1,1,0)
                end if
                
                if (fx(i+1,j,k).ne.covered .and. fy(i+1,j+1,k).ne.covered) then
                else if (fx(i+1,j+1,k).ne.covered .and. fy(i,j+1,k).ne.covered) then
                else
                   flg = clear_neighbor(flg,1,1,0)
                end if
                
                ! x-z
                if (fx(i,j,k).ne.covered .and. fz(i-1,j,k).ne.covered) then
                else if (fx(i,j,k-1).ne.covered .and. fz(i,j,k).ne.covered) then
                else
                   flg = clear_neighbor(flg,-1,0,-1)
                end if
                
                if (fx(i+1,j,k).ne.covered .and. fz(i+1,j,k).ne.covered) then
                else if (fx(i+1,j,k-1).ne.covered .and. fz(i,j,k).ne.covered) then
                else
                   flg = clear_neighbor(flg,1,0,-1)
                end if
                
                if (fx(i,j,k).ne.covered .and. fz(i-1,j,k+1).ne.covered) then
                else if (fx(i,j,k+1).ne.covered .and. fz(i,j,k+1).ne.covered) then
                else
                   flg = clear_neighbor(flg,-1,0,1)
                end if
                
                if (fx(i+1,j,k).ne.covered .and. fz(i+1,j,k+1).ne.covered) then
                else if (fx(i+1,j,k+1).ne.covered .and. fz(i,j,k+1).ne.covered) then
                else
                   flg = clear_neighbor(flg,1,0,1)
                end if
                
                ! y-z
                if (fy(i,j,k).ne.covered .and. fz(i,j-1,k).ne.covered) then
                else if (fy(i,j,k-1).ne.covered .and. fz(i,j,k).ne.covered) then
                else
                   flg = clear_neighbor(flg,0,-1,-1)
                end if
                
                if (fy(i,j+1,k).ne.covered .and. fz(i,j+1,k).ne.covered) then
                else if (fy(i,j+1,k-1).ne.covered .and. fz(i,j,k).ne.covered) then
                else
                   flg = clear_neighbor(flg,0,1,-1)
                end if
                
                if (fy(i,j,k).ne.covered .and. fz(i,j-1,k+1).ne.covered) then
                else if (fy(i,j,k+1).ne.covered .and. fz(i,j,k+1).ne.covered) then
                else
                   flg = clear_neighbor(flg,0,-1,1)
                end if
                
                if (fy(i,j+1,k).ne.covered .and. fz(i,j+1,k+1).ne.covered) then
                else if (fy(i,j+1,k+1).ne.covered .and. fz(i,j,k+1).ne.covered) then
                else
                   flg = clear_neighbor(flg,0,1,1)
                end if

                cell(i,j,k) = flg
             end if
          end do
       end do
    end do

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (.not.is_covered_cell(cell(i,j,k))) then
                flg = cell(i,j,k)
                call get_neighbor_cells(flg, ngbr)
                
                ! -1, -1, -1 corner
                if      (ngbr(-1, 0, 0).eq.1 .and. is_neighbor(cell(i-1,j  ,k  ), 0,-1,-1)) then
                else if (ngbr( 0,-1, 0).eq.1 .and. is_neighbor(cell(i  ,j-1,k  ),-1, 0,-1)) then
                else if (ngbr( 0, 0,-1).eq.1 .and. is_neighbor(cell(i  ,j  ,k-1),-1,-1, 0)) then
                else
                   flg = clear_neighbor(flg, -1,-1,-1)
                end if
                
                ! 1, -1, -1 corner
                if      (ngbr( 1, 0, 0).eq.1 .and. is_neighbor(cell(i+1,j  ,k  ), 0,-1,-1)) then
                else if (ngbr( 0,-1, 0).eq.1 .and. is_neighbor(cell(i  ,j-1,k  ), 1, 0,-1)) then
                else if (ngbr( 0, 0,-1).eq.1 .and. is_neighbor(cell(i  ,j  ,k-1), 1,-1, 0)) then
                else
                   flg = clear_neighbor(flg, 1,-1,-1)
                end if
                
                ! -1, 1, -1 corner
                if      (ngbr(-1, 0, 0).eq.1 .and. is_neighbor(cell(i-1,j  ,k  ), 0, 1,-1)) then
                else if (ngbr( 0, 1, 0).eq.1 .and. is_neighbor(cell(i  ,j+1,k  ),-1, 0,-1)) then
                else if (ngbr( 0, 0,-1).eq.1 .and. is_neighbor(cell(i  ,j  ,k-1),-1, 1, 0)) then
                else
                   flg = clear_neighbor(flg, -1, 1,-1)
                end if
                
                ! 1, 1, -1 corner
                if      (ngbr( 1, 0, 0).eq.1 .and. is_neighbor(cell(i+1,j  ,k  ), 0, 1,-1)) then
                else if (ngbr( 0, 1, 0).eq.1 .and. is_neighbor(cell(i  ,j+1,k  ), 1, 0,-1)) then
                else if (ngbr( 0, 0,-1).eq.1 .and. is_neighbor(cell(i  ,j  ,k-1), 1, 1, 0)) then
                else
                   flg = clear_neighbor(flg, 1, 1,-1)
                end if
                
                ! -1, -1, 1 corner
                if      (ngbr(-1, 0, 0).eq.1 .and. is_neighbor(cell(i-1,j  ,k  ), 0,-1, 1)) then
                else if (ngbr( 0,-1, 0).eq.1 .and. is_neighbor(cell(i  ,j-1,k  ),-1, 0, 1)) then
                else if (ngbr( 0, 0, 1).eq.1 .and. is_neighbor(cell(i  ,j  ,k+1),-1,-1, 0)) then
                else
                   flg = clear_neighbor(flg, -1,-1, 1)
                end if
                
                ! 1, -1, 1 corner
                if      (ngbr( 1, 0, 0).eq.1 .and. is_neighbor(cell(i+1,j  ,k  ), 0,-1, 1)) then
                else if (ngbr( 0,-1, 0).eq.1 .and. is_neighbor(cell(i  ,j-1,k  ), 1, 0, 1)) then
                else if (ngbr( 0, 0, 1).eq.1 .and. is_neighbor(cell(i  ,j  ,k+1), 1,-1, 0)) then
                else
                   flg = clear_neighbor(flg, 1,-1, 1)
                end if
                
                ! -1, 1, 1 corner
                if      (ngbr(-1, 0, 0).eq.1 .and. is_neighbor(cell(i-1,j  ,k  ), 0, 1, 1)) then
                else if (ngbr( 0, 1, 0).eq.1 .and. is_neighbor(cell(i  ,j+1,k  ),-1, 0, 1)) then
                else if (ngbr( 0, 0, 1).eq.1 .and. is_neighbor(cell(i  ,j  ,k+1),-1, 1, 0)) then
                else
                   flg = clear_neighbor(flg, -1,1,1)
                end if
                
                ! 1, 1, 1 corner
                if      (ngbr( 1, 0, 0).eq.1 .and. is_neighbor(cell(i+1,j  ,k  ), 0, 1, 1)) then
                else if (ngbr( 0, 1, 0).eq.1 .and. is_neighbor(cell(i  ,j+1,k  ), 1, 0, 1)) then
                else if (ngbr( 0, 0, 1).eq.1 .and. is_neighbor(cell(i  ,j  ,k+1), 1, 1, 0)) then
                else
                   flg = clear_neighbor(flg, 1,1,1)
                end if
                
                cell(i,j,k) = flg
             end if
          end do
       end do
    end do

  end subroutine amrex_eb2_build_cells
    

  subroutine amrex_eb2_coarsen_from_fine (lo, hi, xlo, xhi, ylo, yhi, zlo, zhi, &
       cvol, cvlo, cvhi, fvol, fvlo, fvhi, ccent, cclo, cchi, fcent, fclo, fchi, &
       cba, cbalo, cbahi, fba, fbalo, fbahi, cbc, cbclo, cbchi, fbc, fbclo, fbchi, &
       cbn, cbnlo, cbnhi, fbn, fbnlo, fbnhi, capx, caxlo, caxhi, fapx, faxlo, faxhi, &
       capy, caylo, cayhi, fapy, faylo, fayhi, capz, cazlo, cazhi, fapz, fazlo, fazhi, &
       cfcx, cfxlo, cfxhi, ffcx, ffxlo, ffxhi, cfcy, cfylo, cfyhi, ffcy, ffylo, ffyhi, &
       cfcz, cfzlo, cfzhi, ffcz, ffzlo, ffzhi, cflag, cflo, cfhi, fflag, fflo, ffhi, ierr) &
       bind(c, name='amrex_eb2_coarsen_from_fine')
    integer, dimension(3), intent(in) :: lo, hi, xlo, xhi, ylo, yhi, zlo, zhi, &
         cvlo, cvhi,  fvlo, fvhi, cclo, cchi, fclo, fchi, &
         cbalo, cbahi, fbalo, fbahi, cbclo, cbchi, fbclo, fbchi, &
         cbnlo, cbnhi, fbnlo, fbnhi, caxlo, caxhi, faxlo, faxhi, &
         caylo, cayhi, faylo, fayhi, cazlo, cazhi, fazlo, fazhi, &
         cfxlo, cfxhi, ffxlo, ffxhi, cfylo, cfyhi, ffylo, ffyhi, &
         cfzlo, cfzhi, ffzlo, ffzhi, cflo, cfhi, fflo, ffhi
    integer, intent(inout) :: ierr
    real(amrex_real), intent(inout) :: cvol ( cvlo(1): cvhi(1), cvlo(2): cvhi(2), cvlo(3): cvhi(3))
    real(amrex_real), intent(in   ) :: fvol ( fvlo(1): fvhi(1), fvlo(2): fvhi(2), fvlo(3): fvhi(3))
    real(amrex_real), intent(inout) :: ccent( cclo(1): cchi(1), cclo(2): cchi(2), cclo(3): cchi(3),3)
    real(amrex_real), intent(in   ) :: fcent( fclo(1): fchi(1), fclo(2): fchi(2), fclo(3): fchi(3),3)
    real(amrex_real), intent(inout) :: cba  (cbalo(1):cbahi(1),cbalo(2):cbahi(2),cbalo(3):cbahi(3))
    real(amrex_real), intent(in   ) :: fba  (fbalo(1):fbahi(1),fbalo(2):fbahi(2),fbalo(3):fbahi(3))
    real(amrex_real), intent(inout) :: cbc  (cbclo(1):cbchi(1),cbclo(2):cbchi(2),cbclo(3):cbchi(3),3)
    real(amrex_real), intent(in   ) :: fbc  (fbclo(1):fbchi(1),fbclo(2):fbchi(2),fbclo(3):fbchi(3),3)
    real(amrex_real), intent(inout) :: cbn  (cbnlo(1):cbnhi(1),cbnlo(2):cbnhi(2),cbnlo(3):cbnhi(3),3)
    real(amrex_real), intent(in   ) :: fbn  (fbnlo(1):fbnhi(1),fbnlo(2):fbnhi(2),fbnlo(3):fbnhi(3),3)
    real(amrex_real), intent(inout) :: capx (caxlo(1):caxhi(1),caxlo(2):caxhi(2),caxlo(3):caxhi(3))
    real(amrex_real), intent(in   ) :: fapx (faxlo(1):faxhi(1),faxlo(2):faxhi(2),faxlo(3):faxhi(3))
    real(amrex_real), intent(inout) :: capy (caylo(1):cayhi(1),caylo(2):cayhi(2),caylo(3):cayhi(3))
    real(amrex_real), intent(in   ) :: fapy (faylo(1):fayhi(1),faylo(2):fayhi(2),faylo(3):fayhi(3))
    real(amrex_real), intent(inout) :: capz (cazlo(1):cazhi(1),cazlo(2):cazhi(2),cazlo(3):cazhi(3))
    real(amrex_real), intent(in   ) :: fapz (fazlo(1):fazhi(1),fazlo(2):fazhi(2),fazlo(3):fazhi(3))
    real(amrex_real), intent(inout) :: cfcx (cfxlo(1):cfxhi(1),cfxlo(2):cfxhi(2),cfxlo(3):cfxhi(3),2)
    real(amrex_real), intent(in   ) :: ffcx (ffxlo(1):ffxhi(1),ffxlo(2):ffxhi(2),ffxlo(3):ffxhi(3),2)
    real(amrex_real), intent(inout) :: cfcy (cfylo(1):cfyhi(1),cfylo(2):cfyhi(2),cfylo(3):cfyhi(3),2)
    real(amrex_real), intent(in   ) :: ffcy (ffylo(1):ffyhi(1),ffylo(2):ffyhi(2),ffylo(3):ffyhi(3),2)
    real(amrex_real), intent(inout) :: cfcz (cfzlo(1):cfzhi(1),cfzlo(2):cfzhi(2),cfzlo(3):cfzhi(3),2)
    real(amrex_real), intent(in   ) :: ffcz (ffzlo(1):ffzhi(1),ffzlo(2):ffzhi(2),ffzlo(3):ffzhi(3),2)
    integer         , intent(inout) :: cflag( cflo(1): cfhi(1), cflo(2): cfhi(2), cflo(3): cfhi(3))
    integer         , intent(in   ) :: fflag( fflo(1): ffhi(1), fflo(2): ffhi(2), fflo(3): ffhi(3))

    integer :: i,j,k, ii,jj,kk, ftype(2,2,2)
    real(amrex_real) :: cvolinv, cbainv, nx, ny, nz, nfac, apinv

    do       k = lo(3), hi(3)
       kk = k*2
       do    j = lo(2), hi(2)
          jj = j*2
          do i = lo(1), hi(1)
             ii = i*2

             ftype = get_cell_type(fflag(ii:ii+1,jj:jj+1,kk:kk+1))
             if (all(ftype.eq.regular_cell)) then
                ! nothing to do
             else if (all(ftype.eq.covered_cell)) then
                call set_covered_cell(cflag(i,j,k))
                cvol(i,j,k) = zero
             else

                call set_single_valued_cell(cflag(i,j,k))

                cvol(i,j,k) = eighth*sum(fvol(ii:ii+1,jj:jj+1,kk:kk+1))
                cvolinv = one/cvol(i,j,k)
                
                ccent(i,j,k,1) = eighth*cvolinv* &
                     ( fvol(ii  ,jj  ,kk  )*(half*fcent(ii  ,jj  ,kk  ,1)-fourth) &
                     + fvol(ii+1,jj  ,kk  )*(half*fcent(ii+1,jj  ,kk  ,1)+fourth) &
                     + fvol(ii  ,jj+1,kk  )*(half*fcent(ii  ,jj+1,kk  ,1)-fourth) &
                     + fvol(ii+1,jj+1,kk  )*(half*fcent(ii+1,jj+1,kk  ,1)+fourth) &
                     + fvol(ii  ,jj  ,kk+1)*(half*fcent(ii  ,jj  ,kk+1,1)-fourth) &
                     + fvol(ii+1,jj  ,kk+1)*(half*fcent(ii+1,jj  ,kk+1,1)+fourth) &
                     + fvol(ii  ,jj+1,kk+1)*(half*fcent(ii  ,jj+1,kk+1,1)-fourth) &
                     + fvol(ii+1,jj+1,kk+1)*(half*fcent(ii+1,jj+1,kk+1,1)+fourth) )
                ccent(i,j,k,2) = eighth*cvolinv* &
                     ( fvol(ii  ,jj  ,kk  )*(half*fcent(ii  ,jj  ,kk  ,2)-fourth) &
                     + fvol(ii+1,jj  ,kk  )*(half*fcent(ii+1,jj  ,kk  ,2)-fourth) &
                     + fvol(ii  ,jj+1,kk  )*(half*fcent(ii  ,jj+1,kk  ,2)+fourth) &
                     + fvol(ii+1,jj+1,kk  )*(half*fcent(ii+1,jj+1,kk  ,2)+fourth) &
                     + fvol(ii  ,jj  ,kk+1)*(half*fcent(ii  ,jj  ,kk+1,2)-fourth) &
                     + fvol(ii+1,jj  ,kk+1)*(half*fcent(ii+1,jj  ,kk+1,2)-fourth) &
                     + fvol(ii  ,jj+1,kk+1)*(half*fcent(ii  ,jj+1,kk+1,2)+fourth) &
                     + fvol(ii+1,jj+1,kk+1)*(half*fcent(ii+1,jj+1,kk+1,2)+fourth) )
                ccent(i,j,k,3) = eighth*cvolinv* &
                     ( fvol(ii  ,jj  ,kk  )*(half*fcent(ii  ,jj  ,kk  ,3)-fourth) &
                     + fvol(ii+1,jj  ,kk  )*(half*fcent(ii+1,jj  ,kk  ,3)-fourth) &
                     + fvol(ii  ,jj+1,kk  )*(half*fcent(ii  ,jj+1,kk  ,3)-fourth) &
                     + fvol(ii+1,jj+1,kk  )*(half*fcent(ii+1,jj+1,kk  ,3)-fourth) &
                     + fvol(ii  ,jj  ,kk+1)*(half*fcent(ii  ,jj  ,kk+1,3)+fourth) &
                     + fvol(ii+1,jj  ,kk+1)*(half*fcent(ii+1,jj  ,kk+1,3)+fourth) &
                     + fvol(ii  ,jj+1,kk+1)*(half*fcent(ii  ,jj+1,kk+1,3)+fourth) &
                     + fvol(ii+1,jj+1,kk+1)*(half*fcent(ii+1,jj+1,kk+1,3)+fourth) )
                
                cba(i,j,k) = fourth*sum(fba(ii:ii+1,jj:jj+1,kk:kk+1))
                cbainv = one/cba(i,j,k)
                
                cbc(i,j,k,1) = fourth*cbainv* &
                     ( fba(ii  ,jj  ,kk  )*(half*fbc(ii  ,jj  ,kk  ,1)-fourth) &
                     + fba(ii+1,jj  ,kk  )*(half*fbc(ii+1,jj  ,kk  ,1)+fourth) &
                     + fba(ii  ,jj+1,kk  )*(half*fbc(ii  ,jj+1,kk  ,1)-fourth) &
                     + fba(ii+1,jj+1,kk  )*(half*fbc(ii+1,jj+1,kk  ,1)+fourth) &
                     + fba(ii  ,jj  ,kk+1)*(half*fbc(ii  ,jj  ,kk+1,1)-fourth) &
                     + fba(ii+1,jj  ,kk+1)*(half*fbc(ii+1,jj  ,kk+1,1)+fourth) &
                     + fba(ii  ,jj+1,kk+1)*(half*fbc(ii  ,jj+1,kk+1,1)-fourth) &
                     + fba(ii+1,jj+1,kk+1)*(half*fbc(ii+1,jj+1,kk+1,1)+fourth) )
                cbc(i,j,k,2) = fourth*cbainv* &
                     ( fba(ii  ,jj  ,kk  )*(half*fbc(ii  ,jj  ,kk  ,2)-fourth) &
                     + fba(ii+1,jj  ,kk  )*(half*fbc(ii+1,jj  ,kk  ,2)-fourth) &
                     + fba(ii  ,jj+1,kk  )*(half*fbc(ii  ,jj+1,kk  ,2)+fourth) &
                     + fba(ii+1,jj+1,kk  )*(half*fbc(ii+1,jj+1,kk  ,2)+fourth) &
                     + fba(ii  ,jj  ,kk+1)*(half*fbc(ii  ,jj  ,kk+1,2)-fourth) &
                     + fba(ii+1,jj  ,kk+1)*(half*fbc(ii+1,jj  ,kk+1,2)-fourth) &
                     + fba(ii  ,jj+1,kk+1)*(half*fbc(ii  ,jj+1,kk+1,2)+fourth) &
                     + fba(ii+1,jj+1,kk+1)*(half*fbc(ii+1,jj+1,kk+1,2)+fourth) )
                cbc(i,j,k,3) = fourth*cbainv* &
                     ( fba(ii  ,jj  ,kk  )*(half*fbc(ii  ,jj  ,kk  ,3)-fourth) &
                     + fba(ii+1,jj  ,kk  )*(half*fbc(ii+1,jj  ,kk  ,3)-fourth) &
                     + fba(ii  ,jj+1,kk  )*(half*fbc(ii  ,jj+1,kk  ,3)-fourth) &
                     + fba(ii+1,jj+1,kk  )*(half*fbc(ii+1,jj+1,kk  ,3)-fourth) &
                     + fba(ii  ,jj  ,kk+1)*(half*fbc(ii  ,jj  ,kk+1,3)+fourth) &
                     + fba(ii+1,jj  ,kk+1)*(half*fbc(ii+1,jj  ,kk+1,3)+fourth) &
                     + fba(ii  ,jj+1,kk+1)*(half*fbc(ii  ,jj+1,kk+1,3)+fourth) &
                     + fba(ii+1,jj+1,kk+1)*(half*fbc(ii+1,jj+1,kk+1,3)+fourth) )
                
                nx =   fbn(ii  ,jj  ,kk  ,1)*fba(ii  ,jj  ,kk  ) &
                     + fbn(ii+1,jj  ,kk  ,1)*fba(ii+1,jj  ,kk  ) &
                     + fbn(ii  ,jj+1,kk  ,1)*fba(ii  ,jj+1,kk  ) &
                     + fbn(ii+1,jj+1,kk  ,1)*fba(ii+1,jj+1,kk  ) &
                     + fbn(ii  ,jj  ,kk+1,1)*fba(ii  ,jj  ,kk+1) &
                     + fbn(ii+1,jj  ,kk+1,1)*fba(ii+1,jj  ,kk+1) &
                     + fbn(ii  ,jj+1,kk+1,1)*fba(ii  ,jj+1,kk+1) &
                     + fbn(ii+1,jj+1,kk+1,1)*fba(ii+1,jj+1,kk+1)
                ny =   fbn(ii  ,jj  ,kk  ,2)*fba(ii  ,jj  ,kk  ) &
                     + fbn(ii+1,jj  ,kk  ,2)*fba(ii+1,jj  ,kk  ) &
                     + fbn(ii  ,jj+1,kk  ,2)*fba(ii  ,jj+1,kk  ) &
                     + fbn(ii+1,jj+1,kk  ,2)*fba(ii+1,jj+1,kk  ) &
                     + fbn(ii  ,jj  ,kk+1,2)*fba(ii  ,jj  ,kk+1) &
                     + fbn(ii+1,jj  ,kk+1,2)*fba(ii+1,jj  ,kk+1) &
                     + fbn(ii  ,jj+1,kk+1,2)*fba(ii  ,jj+1,kk+1) &
                     + fbn(ii+1,jj+1,kk+1,2)*fba(ii+1,jj+1,kk+1)
                nz =   fbn(ii  ,jj  ,kk  ,3)*fba(ii  ,jj  ,kk  ) &
                     + fbn(ii+1,jj  ,kk  ,3)*fba(ii+1,jj  ,kk  ) &
                     + fbn(ii  ,jj+1,kk  ,3)*fba(ii  ,jj+1,kk  ) &
                     + fbn(ii+1,jj+1,kk  ,3)*fba(ii+1,jj+1,kk  ) &
                     + fbn(ii  ,jj  ,kk+1,3)*fba(ii  ,jj  ,kk+1) &
                     + fbn(ii+1,jj  ,kk+1,3)*fba(ii+1,jj  ,kk+1) &
                     + fbn(ii  ,jj+1,kk+1,3)*fba(ii  ,jj+1,kk+1) &
                     + fbn(ii+1,jj+1,kk+1,3)*fba(ii+1,jj+1,kk+1)
                if (nx.eq.zero .and. ny.eq.zero .and. nz.eq.zero) then
                   ierr = 1
                   return
                end if
                nfac = one/sqrt(nx*nx+ny*ny+nz*nz)
                cbn(i,j,k,1) = nx*nfac
                cbn(i,j,k,2) = ny*nfac
                cbn(i,j,k,3) = nz*nfac
             end if
          end do
       end do
    end do

    do       k = xlo(3), xhi(3)
       kk = k*2
       do    j = xlo(2), xhi(2)
          jj = j*2
          do i = xlo(1), xhi(1)
             ii = i*2

             capx(i,j,k) = fourth*(fapx(ii,jj,kk)+fapx(ii,jj+1,kk)+fapx(ii,jj,kk+1)+fapx(ii,jj+1,kk+1))
             if (capx(i,j,k) .ne. zero) then
                apinv = one/capx(i,j,k)
                cfcx(i,j,k,1) = fourth*apinv* &
                     ( fapx(ii,jj  ,kk  )*(half*ffcx(ii,jj  ,kk  ,1)-fourth) &
                     + fapx(ii,jj+1,kk  )*(half*ffcx(ii,jj+1,kk  ,1)+fourth) &
                     + fapx(ii,jj  ,kk+1)*(half*ffcx(ii,jj  ,kk+1,1)-fourth) &
                     + fapx(ii,jj+1,kk+1)*(half*ffcx(ii,jj+1,kk+1,1)+fourth) )
                cfcx(i,j,k,2) = fourth*apinv* &
                     ( fapx(ii,jj  ,kk  )*(half*ffcx(ii,jj  ,kk  ,2)-fourth) &
                     + fapx(ii,jj+1,kk  )*(half*ffcx(ii,jj+1,kk  ,2)-fourth) &
                     + fapx(ii,jj  ,kk+1)*(half*ffcx(ii,jj  ,kk+1,2)+fourth) &
                     + fapx(ii,jj+1,kk+1)*(half*ffcx(ii,jj+1,kk+1,2)+fourth) )                
             end if
          end do
       end do
    end do

    do       k = ylo(3), yhi(3)
       kk = k*2
       do    j = ylo(2), yhi(2)
          jj = j*2
          do i = ylo(1), yhi(1)
             ii = i*2

             capy(i,j,k) = fourth*(fapy(ii,jj,kk)+fapy(ii+1,jj,kk)+fapy(ii,jj,kk+1)+fapy(ii+1,jj,kk+1))
             if (capy(i,j,k) .ne. zero) then
                apinv = one/capy(i,j,k)
                cfcy(i,j,k,1) = fourth*apinv* &
                     ( fapy(ii  ,jj,kk  )*(half*ffcy(ii  ,jj,kk  ,1)-fourth) &
                     + fapy(ii+1,jj,kk  )*(half*ffcy(ii+1,jj,kk  ,1)+fourth) &
                     + fapy(ii  ,jj,kk+1)*(half*ffcy(ii  ,jj,kk+1,1)-fourth) &
                     + fapy(ii+1,jj,kk+1)*(half*ffcy(ii+1,jj,kk+1,1)+fourth) )
                cfcy(i,j,k,2) = fourth*apinv* &
                     ( fapy(ii  ,jj,kk  )*(half*ffcy(ii  ,jj,kk  ,2)-fourth) &
                     + fapy(ii+1,jj,kk  )*(half*ffcy(ii+1,jj,kk  ,2)-fourth) &
                     + fapy(ii  ,jj,kk+1)*(half*ffcy(ii  ,jj,kk+1,2)+fourth) &
                     + fapy(ii+1,jj,kk+1)*(half*ffcy(ii+1,jj,kk+1,2)+fourth) )                
             end if
          end do
       end do
    end do

    do       k = zlo(3), zhi(3)
       kk = k*2
       do    j = zlo(2), zhi(2)
          jj = j*2
          do i = zlo(1), zhi(1)
             ii = i*2

             capz(i,j,k) = fourth*(fapz(ii,jj,kk)+fapz(ii+1,jj,kk)+fapz(ii,jj+1,kk)+fapz(ii+1,jj+1,kk))
             if (capz(i,j,k) .ne. zero) then
                apinv = one/capz(i,j,k)
                cfcz(i,j,k,1) = fourth*apinv* &
                     ( fapz(ii  ,jj  ,kk)*(half*ffcz(ii  ,jj  ,kk,1)-fourth) &
                     + fapz(ii+1,jj  ,kk)*(half*ffcz(ii+1,jj  ,kk,1)+fourth) &
                     + fapz(ii  ,jj+1,kk)*(half*ffcz(ii  ,jj+1,kk,1)-fourth) &
                     + fapz(ii+1,jj+1,kk)*(half*ffcz(ii+1,jj+1,kk,1)+fourth) )
                cfcz(i,j,k,2) = fourth*apinv* &
                     ( fapz(ii  ,jj  ,kk)*(half*ffcz(ii  ,jj  ,kk,2)-fourth) &
                     + fapz(ii+1,jj  ,kk)*(half*ffcz(ii+1,jj  ,kk,2)-fourth) &
                     + fapz(ii  ,jj+1,kk)*(half*ffcz(ii  ,jj+1,kk,2)+fourth) &
                     + fapz(ii+1,jj+1,kk)*(half*ffcz(ii+1,jj+1,kk,2)+fourth) )                
             end if
          end do
       end do
    end do

  end subroutine amrex_eb2_coarsen_from_fine


  subroutine amrex_eb2_build_cellflag_from_ap (lo, hi, &
       apx, xlo, xhi, apy, ylo, yhi, apz, zlo, zhi, cflag, flo, fhi) &
       bind(c,name='amrex_eb2_build_cellflag_from_ap')
    integer, dimension(3) :: lo, hi, xlo, xhi, ylo, yhi, zlo, zhi, flo, fhi
    real(amrex_real), intent(in   ) ::  apx (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) ::  apy (ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    real(amrex_real), intent(in   ) ::  apz (zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))
    integer         , intent(inout) :: cflag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

    integer :: i,j,k, flg

    ! By default, all neighbors are already set.
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             cflag(i,j,k) = clear_allneighbors(cflag(i,j,k))
                
             if (.not.is_covered_cell(cflag(i,j,k))) then
                flg = clear_allneighbors(cflag(i,j,k))
                
                if (apx(i  ,j,k).ne.zero) flg = set_neighbor(flg, -1,  0,  0)
                if (apx(i+1,j,k).ne.zero) flg = set_neighbor(flg,  1,  0,  0)
                if (apy(i,j  ,k).ne.zero) flg = set_neighbor(flg,  0, -1,  0)
                if (apy(i,j+1,k).ne.zero) flg = set_neighbor(flg,  0,  1,  0)
                if (apz(i,j,k  ).ne.zero) flg = set_neighbor(flg,  0,  0, -1)
                if (apz(i,j,k+1).ne.zero) flg = set_neighbor(flg,  0,  0,  1)
                
                if ( (apx(i,j,k).ne.zero .and. apy(i-1,j,k).ne.zero) .or. &
                     (apy(i,j,k).ne.zero .and. apx(i,j-1,k).ne.zero) ) then
                   flg = set_neighbor(flg, -1, -1, 0)
                   if (apz(i-1,j-1,k  ).ne.zero) flg = set_neighbor(flg,-1,-1,-1)
                   if (apz(i-1,j-1,k+1).ne.zero) flg = set_neighbor(flg,-1,-1, 1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apy(i+1,j,k).ne.zero) .or. &
                     (apy(i,j,k).ne.zero .and. apx(i+1,j-1,k).ne.zero) ) then
                   flg = set_neighbor(flg, 1, -1, 0)
                   if (apz(i+1,j-1,k  ).ne.zero) flg = set_neighbor(flg, 1,-1,-1)
                   if (apz(i+1,j-1,k+1).ne.zero) flg = set_neighbor(flg, 1,-1, 1)
                end if
                
                if ( (apx(i,j,k).ne.zero .and. apy(i-1,j+1,k).ne.zero) .or. &
                     (apy(i,j+1,k).ne.zero .and. apx(i,j+1,k).ne.zero) ) then
                   flg = set_neighbor(flg, -1, 1, 0)
                   if (apz(i-1,j+1,k  ).ne.zero) flg = set_neighbor(flg,-1, 1,-1)
                   if (apz(i-1,j+1,k+1).ne.zero) flg = set_neighbor(flg,-1, 1, 1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apy(i+1,j+1,k).ne.zero) .or. &
                     (apy(i,j+1,k).ne.zero .and. apx(i+1,j+1,k).ne.zero) ) then
                   flg = set_neighbor(flg, 1, 1, 0)
                   if (apz(i+1,j+1,k  ).ne.zero) flg = set_neighbor(flg, 1, 1,-1)
                   if (apz(i+1,j+1,k+1).ne.zero) flg = set_neighbor(flg, 1, 1, 1)
                end if
                
                if ( (apx(i,j,k).ne.zero .and. apz(i-1,j,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apx(i,j,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, -1, 0, -1)
                   if (apy(i-1,j  ,k-1).ne.zero) flg = set_neighbor(flg,-1,-1,-1)
                   if (apy(i-1,j+1,k-1).ne.zero) flg = set_neighbor(flg,-1, 1,-1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apz(i+1,j,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apx(i+1,j,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, 1, 0, -1)
                   if (apy(i+1,j  ,k-1).ne.zero) flg = set_neighbor(flg, 1,-1,-1)
                   if (apy(i+1,j+1,k-1).ne.zero) flg = set_neighbor(flg, 1, 1,-1)
                end if
                
                if ( (apx(i,j,k).ne.zero .and. apz(i-1,j,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apx(i,j,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, -1, 0, 1)
                   if (apy(i-1,j  ,k+1).ne.zero) flg = set_neighbor(flg,-1,-1, 1)
                   if (apy(i-1,j+1,k+1).ne.zero) flg = set_neighbor(flg,-1, 1, 1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apz(i+1,j,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apx(i+1,j,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, 1, 0, 1)
                   if (apy(i+1,j  ,k+1).ne.zero) flg = set_neighbor(flg, 1,-1, 1)
                   if (apy(i+1,j+1,k+1).ne.zero) flg = set_neighbor(flg, 1, 1, 1)
                end if
                
                if ( (apy(i,j,k).ne.zero .and. apz(i,j-1,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apy(i,j,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, -1, -1)
                   if (apx(i  ,j-1,k-1).ne.zero) flg = set_neighbor(flg,-1,-1,-1)
                   if (apx(i+1,j-1,k-1).ne.zero) flg = set_neighbor(flg, 1,-1,-1)
                end if
                
                if ( (apy(i,j+1,k).ne.zero .and. apz(i,j+1,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apy(i,j+1,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, 1, -1)
                   if (apx(i  ,j+1,k-1).ne.zero) flg = set_neighbor(flg,-1, 1,-1)
                   if (apx(i+1,j+1,k-1).ne.zero) flg = set_neighbor(flg, 1, 1,-1)
                end if
                
                if ( (apy(i,j,k).ne.zero .and. apz(i,j-1,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apy(i,j,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, -1, 1)
                   if (apx(i  ,j-1,k+1).ne.zero) flg = set_neighbor(flg,-1,-1, 1)
                   if (apx(i+1,j-1,k+1).ne.zero) flg = set_neighbor(flg, 1,-1, 1)
                end if
                
                if ( (apy(i,j+1,k).ne.zero .and. apz(i,j+1,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apy(i,j+1,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, 1, 1)
                   if (apx(i  ,j+1,k+1).ne.zero) flg = set_neighbor(flg,-1, 1, 1)
                   if (apx(i+1,j+1,k+1).ne.zero) flg = set_neighbor(flg, 1, 1, 1)
                end if
                   
                cflag(i,j,k) = flg
             end if
          end do
       end do
    end do
  end subroutine amrex_eb2_build_cellflag_from_ap


  subroutine amrex_eb2_check_mvmc (cclo, cchi, ndlo, ndhi, cls, clo, chi, fls, flo, fhi, &
       ncuts, tlo, thi, ierr) &
       bind(c,name='amrex_eb2_check_mvmc')
    integer, dimension(3), intent(in) :: cclo,cchi,ndlo,ndhi,clo,chi,flo,fhi,tlo,thi
    real(amrex_real), intent(inout) :: cls(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(amrex_real), intent(in   ) :: fls(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    integer, intent(inout) :: ncuts(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),6)
    integer, intent(inout) :: ierr

    integer :: i,j,k, ii,jj,kk, n, nopen

    do       k = ndlo(3), ndhi(3)
       kk = k*2
       do    j = ndlo(2), ndhi(2)
          jj = j*2
          do i = ndlo(1), ndhi(1)
             ii = i*2
             cls(i,j,k) = fls(ii,jj,kk)
          end do
       end do
    end do

    ! x-edges
    do       k = cclo(3), cchi(3)+1
       kk = k*2
       do    j = cclo(2), cchi(2)+1
          jj = j*2
          do i = cclo(1), cchi(1)
             ii = i*2
             ncuts(i,j,k,1) = 0
             if (has_cut(fls(ii,jj,kk),fls(ii+1,jj,kk))) then
                ncuts(i,j,k,1) = ncuts(i,j,k,1)+1
             end if
             if (has_cut(fls(ii+1,jj,kk),fls(ii+2,jj,kk))) then
                ncuts(i,j,k,1) = ncuts(i,j,k,1)+1
             end if
             if (ncuts(i,j,k,1) .eq. 2) then
                ierr = 1
                return
             end if
          end do
       end do
    end do

    ! y-edges
    do       k = cclo(3), cchi(3)+1
       kk = k*2
       do    j = cclo(2), cchi(2)
          jj = j*2
          do i = cclo(1), cchi(1)+1
             ii = i*2
             ncuts(i,j,k,2) = 0
             if (has_cut(fls(ii,jj,kk),fls(ii,jj+1,kk))) then
                ncuts(i,j,k,2) = ncuts(i,j,k,2)+1
             end if
             if (has_cut(fls(ii,jj+1,kk),fls(ii,jj+2,kk))) then
                ncuts(i,j,k,2) = ncuts(i,j,k,2)+1
             end if
             if (ncuts(i,j,k,2) .eq. 2) then
                ierr = 1
                return
             end if
          end do
       end do
    end do

    ! z-edges
    do       k = cclo(3), cchi(3)
       kk = k*2
       do    j = cclo(2), cchi(2)+1
          jj = j*2
          do i = cclo(1), cchi(1)+1
             ii = i*2
             ncuts(i,j,k,3) = 0
             if (has_cut(fls(ii,jj,kk),fls(ii,jj,kk+1))) then
                ncuts(i,j,k,3) = ncuts(i,j,k,3)+1
             end if
             if (has_cut(fls(ii,jj,kk+1),fls(ii,jj,kk+2))) then
                ncuts(i,j,k,3) = ncuts(i,j,k,3)+1
             end if
             if (ncuts(i,j,k,3) .eq. 2) then
                ierr = 1
                return
             end if
          end do
       end do
    end do

    ! x-faces
    do       k = cclo(3), cchi(3)
       do    j = cclo(2), cchi(2)
          do i = cclo(1), cchi(1)+1
             n =    ncuts(i,j,k,2) + ncuts(i,j,k+1,2) &
                  + ncuts(i,j,k,3) + ncuts(i,j+1,k,3)
             if (n .eq. 0) then
                ncuts(i,j,k,4) = 0
             else if (n .eq. 2) then
                ncuts(i,j,k,4) = 1
             else if (n .eq. 4) then
                ierr = 1
                return
             else
                call amrex_error("amrex_eb2_check_mvmc: how did this happen? wrong nubmer of cuts on x-face", n)
             end if
          end do
       end do
    end do
    
    ! y-faces
    do       k = cclo(3), cchi(3)
       do    j = cclo(2), cchi(2)+1
          do i = cclo(1), cchi(1)
             n =    ncuts(i,j,k,1) + ncuts(i,j,k+1,1) &
                  + ncuts(i,j,k,3) + ncuts(i+1,j,k,3)
             if (n .eq. 0) then
                ncuts(i,j,k,5) = 0
             else if (n .eq. 2) then
                ncuts(i,j,k,5) = 1
             else if (n .eq. 4) then
                ierr = 1
                return
             else
                call amrex_error("amrex_eb2_check_mvmc: how did this happen? wrong nubmer of cuts on y-face", n)
             end if
          end do
       end do
    end do

    ! z-faces
    do       k = cclo(3), cchi(3)+1
       do    j = cclo(2), cchi(2)
          do i = cclo(1), cchi(1)
             n =    ncuts(i,j,k,1) + ncuts(i,j+1,k,1) &
                  + ncuts(i,j,k,2) + ncuts(i+1,j,k,2)
             if (n .eq. 0) then
                ncuts(i,j,k,6) = 0
             else if (n .eq. 2) then
                ncuts(i,j,k,6) = 1
             else if (n .eq. 4) then
                ierr = 1
                return
             else
                call amrex_error("amrex_eb2_check_mvmc: how did this happen? wrong nubmer of cuts on x-face", n)
             end if
          end do
       end do
    end do

    do       k = cclo(3), cchi(3)
       do    j = cclo(2), cchi(2)
          do i = cclo(1), cchi(1)
             if ( ncuts(i  ,j  ,k  ,4).eq.1 .and. &
                  ncuts(i+1,j  ,k  ,4).eq.1 .and. &
                  ncuts(i  ,j  ,k  ,5).eq.1 .and. &
                  ncuts(i  ,j+1,k  ,5).eq.1 .and. &
                  ncuts(i  ,j  ,k  ,6).eq.1 .and. &
                  ncuts(i  ,j  ,k+1,6).eq.1 ) then
                ii = i*2
                jj = j*2
                kk = k*2
                nopen = 0
                if (fls(ii  ,jj  ,kk  ) .lt. zero) nopen = nopen + 1
                if (fls(ii+2,jj  ,kk  ) .lt. zero) nopen = nopen + 1
                if (fls(ii  ,jj+2,kk  ) .lt. zero) nopen = nopen + 1
                if (fls(ii+2,jj+2,kk  ) .lt. zero) nopen = nopen + 1
                if (fls(ii  ,jj  ,kk+2) .lt. zero) nopen = nopen + 1
                if (fls(ii+2,jj  ,kk+2) .lt. zero) nopen = nopen + 1
                if (fls(ii  ,jj+2,kk+2) .lt. zero) nopen = nopen + 1
                if (fls(ii+2,jj+2,kk+2) .lt. zero) nopen = nopen + 1
                if (nopen .eq. 2 .or. nopen .eq. 6) then
                   ierr = 1
                   return
                else if (nopen .ne. 4) then
                   call amrex_error("amrex_eb2_check_mvmc: how did this happen? nopen", nopen)
                end if
             end if
          end do
       end do
    end do
    
  contains
    pure function has_cut(a,b)
      real(amrex_real), intent(in) :: a,b
      logical has_cut
      has_cut = (a.ge.zero .and. b.lt.zero) .or. (b.ge.zero .and. a.lt.zero)
    end function has_cut

  end subroutine amrex_eb2_check_mvmc
  
end module amrex_eb2_3d_module
