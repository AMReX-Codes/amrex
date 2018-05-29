
module amrex_eb2_3d_module

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module, only : zero, one, two, three, four, six, seven, eight, nine, ten, &
       half, third, fourth, sixth, eighth, twelfth
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
                cell(i,j,k) = covered
             else if (  s(i,j  ,k  ).lt.zero .and. s(i+1,j  ,k  ).lt.zero &
                  .and. s(i,j+1,k  ).lt.zero .and. s(i+1,j+1,k  ).lt.zero &
                  .and. s(i,j  ,k+1).lt.zero .and. s(i+1,j  ,k+1).lt.zero &
                  .and. s(i,j+1,k+1).lt.zero .and. s(i+1,j+1,k+1).lt.zero ) then
                cell(i,j,k) = regular
             else
                cell(i,j,k) = irregular
             end if
          end do
       end do
    end do

    ! x-face
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+2
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
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+2
          do i = lo(1)-1, hi(1)+1
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
    do       k = lo(3)-1, hi(3)+2
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
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
    do       k = lo(3)-1, hi(3)+2
       do    j = lo(2)-1, hi(2)+2
          do i = lo(1)-1, hi(1)+1
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
    do       k = lo(3)-1, hi(3)+2
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+2
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
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+2
          do i = lo(1)-1, hi(1)+2
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
    real(amrex_real), intent(inout) ::    fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),3)
    real(amrex_real), intent(inout) ::    fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),3)
    real(amrex_real), intent(inout) ::    fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),3)
    real(amrex_real), intent(inout) ::    m2x(mxlo(1):mxhi(1),mxlo(2):mxhi(2),mxlo(3):mxhi(3),3)
    real(amrex_real), intent(inout) ::    m2y(mylo(1):myhi(1),mylo(2):myhi(2),mylo(3):myhi(3),3)
    real(amrex_real), intent(inout) ::    m2z(mzlo(1):mzhi(1),mzlo(2):mzhi(2),mzlo(3):mzhi(3),3)

    integer :: i,j,k, ncuts
    real(amrex_real) :: cut
    real(amrex_real) :: lxm, lxp, lym, lyp, lzm, lzp, bcx, bcy, bcz
    real(amrex_real), parameter :: twentyfourth = 1.d0/24.d0

    ! x-face
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+2
             if (fx(i,j,k) .eq. regular) then
                apx(i,j,k) = one
                fcx(i,j,k,:) = zero
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
                   call amrex_error("amrex_eb2_build_faces: more than 2 cuts not suported")
                else if (ncuts .lt. 2) then
                   print *, "ncuts = ", i,j,k,ncuts
                   call amrex_error("amrex_eb2_build_faces: irregular face with less than 2 cuts???")
                end if

                if (lym.eq.lyp .and. lzm.eq.lzp) then
                   apx(i,j,k) = one
                   fcx(i,j,k,:) = zero
                   m2x(i,j,k,1) = twelfth
                   m2x(i,j,k,2) = twelfth
                   m2x(i,j,k,3) = zero
                else
                   fcx(i,j,k,1) = zero
                   bcy = half*bcy - half
                   bcz = half*bcz - half
                   call cut_face_2d(apx(i,j,k),fcx(i,j,k,2),fcx(i,j,k,3), &
                        m2x(i,j,k,1),m2x(i,j,k,2),m2x(i,j,k,3),lzm,lzp,lym,lyp,bcy,bcz)
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
                apy(i,j,k) = one
                fcy(i,j,k,:) = zero
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
                   call amrex_error("amrex_eb2_build_faces: more than 2 cuts not supported")
                else if (ncuts .lt. 2) then
                   print *, "ncuts = ", i,j,k,ncuts
                   call amrex_error("amrex_eb2_build_faces: irregular face with less than 2 cuts???")
                end if

                if (lxm.eq.lxp .and. lzm.eq.lzp) then
                   apy(i,j,k) = one
                   fcy(i,j,k,:) = zero
                   m2y(i,j,k,1) = twelfth
                   m2y(i,j,k,2) = twelfth
                   m2y(i,j,k,3) = zero
                else
                   fcy(i,j,k,2) = zero
                   bcx = half*bcx - half
                   bcz = half*bcz - half
                   call cut_face_2d(apy(i,j,k),fcy(i,j,k,1),fcy(i,j,k,3),&
                        m2y(i,j,k,1),m2y(i,j,k,2),m2y(i,j,k,3),lzm,lzp,lxm,lxp,bcx,bcz)
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
                apz(i,j,k) = one
                fcz(i,j,k,:) = zero
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
                   call amrex_error("amrex_eb2_build_faces: more than 2 cuts not supported")
                else if (ncuts .lt. 2) then
                   print *, "ncuts = ", i,j,k,ncuts
                   call amrex_error("amrex_eb2_build_faces: irregular face with less than 2 cuts???")
                end if

                if (lxm.eq.lxp .and. lym.eq.lyp) then
                   apz(i,j,k) = one
                   fcz(i,j,k,:) = zero
                   m2z(i,j,k,1) = twelfth
                   m2z(i,j,k,2) = twelfth
                   m2z(i,j,k,3) = zero
                else
                   fcz(i,j,k,3) = zero
                   bcx = half*bcx - half
                   bcy = half*bcy - half
                   call cut_face_2d(apz(i,j,k),fcz(i,j,k,1),fcz(i,j,k,2), &
                        m2z(i,j,k,1),m2z(i,j,k,2),m2z(i,j,k,3),lym,lyp,lxm,lxp,bcx,bcy)
                end if
             end if
          end do
       end do
    end do

    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             if (cell(i,j,k) .eq. irregular) then
                if  (apx(i,j,k).eq.zero .and. apx(i+1,j,k).eq.zero .and. &
                     apy(i,j,k).eq.zero .and. apy(i,j+1,k).eq.zero .and. &
                     apz(i,j,k).eq.zero .and. apz(i,j,k+1).eq.zero) then
                   cell(i,j,k) = covered
                   fx(i,j,k) = covered
                   fx(i+1,j,k) = covered
                   fy(i,j,k) = covered
                   fy(i,j+1,k) = covered
                   fz(i,j,k) = covered
                   fz(i,j,k+1) = covered
                else if (apx(i,j,k).eq.one .and. apx(i+1,j,k).eq.one .and. &
                     &   apy(i,j,k).eq.one .and. apy(i,j+1,k).eq.one .and. &
                     &   apz(i,j,k).eq.one .and. apz(i,j,k+1).eq.one) then
                   cell(i,j,k) = regular
                   fx(i,j,k) = regular
                   fx(i+1,j,k) = regular
                   fy(i,j,k) = regular
                   fy(i,j+1,k) = regular
                   fz(i,j,k) = regular
                   fz(i,j,k+1) = regular
                end if
             end if
          end do
       end do
    end do

  contains

    subroutine cut_face_2d (areafrac, centx, centy, Sx2, Sy2, Sxy, axm, axp, aym, ayp, bcx, bcy)
      real(amrex_real), intent(out) :: areafrac, centx, centy, Sx2, Sy2, Sxy
      real(amrex_real), intent(in) :: axm, axp, aym, ayp, bcx, bcy

      real(amrex_real) :: apnorm, apnorminv, anrmx, anrmy, x_ym, x_yp, y_xm, y_xp, aa, af1, af2
      real(amrex_real) :: rhs(8), b(7), signx, signy, den, ny2, ny3, ny4, ny5

      apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2)
      apnorminv = one/apnorm
      anrmx = (axm-axp) * apnorminv  ! pointing to the wall
      anrmy = (aym-ayp) * apnorminv

      if (anrmx .eq. zero) then
         areafrac = axm
         centx = zero
         centy = (eighth*(ayp-aym) + anrmy*half*bcy**2)/areafrac
         Sx2 = twentyfourth*(axm+axp)
         Sy2 = twentyfourth*(ayp+aym) + anrmy*third*bcy**3
         Sxy = zero
      else if (anrmy .eq. zero) then
         areafrac = aym
         centx = (eighth*(axp-axm) + anrmx*half*bcx**2)/areafrac
         centy = zero
         Sx2 = twentyfourth*(axp+axm) + anrmx*third*bcx**3
         Sy2 = twentyfourth*(ayp+aym)
         Sxy = zero
      else
         if (anrmx .gt. zero) then
            x_ym = -half + aym
            x_yp = -half + ayp
            signx = one
         else
            x_ym = half - aym
            x_yp = half - ayp
            signx = -one
         end if
         aa = abs(anrmx)/anrmy
         af1 = half*(axm+axp) + aa*half*(x_ym**2-x_yp**2)
         centx = eighth*(axp-axm) + aa*sixth*(x_ym**3-x_yp**3)

         if (anrmy .gt. zero) then
            y_xm = -half + axm
            y_xp = -half + axp
            signy = one
         else
            y_xm = half - axm
            y_xp = half - axp
            signy = -one
         end if
         aa = abs(anrmy)/anrmx
         af2 = half*(aym+ayp) + aa*half*(y_xm**2-y_xp**2)
         centy = eighth*(ayp-aym) + aa*sixth*(y_xm**3-y_xp**3)

         rhs(1) = signx*fourth*(x_ym**4-x_yp**4)
         rhs(2) = -twentyfourth - signx*sixth*(x_ym**3+x_yp**3)
         rhs(3) = signx*eighth*(x_ym**2-x_yp**2)
         rhs(4) = -eighth*(aym+ayp)
         rhs(5) = eighth*(axm+axp)
         rhs(6) = -signx*eighth*(y_xm**2-y_xp**2)
         rhs(7) = twentyfourth + signx*sixth*(y_xm**3+y_xp**3)
         rhs(8) = -signx*fourth*(y_xm**4-y_xp**4)

         b(1) = -rhs(2) + three*rhs(5)
         b(2) = -three*rhs(4) + rhs(7)
         b(3) = two*(-rhs(3) + rhs(6))
         b(4) = anrmy*rhs(1) - anrmx*rhs(5)
         b(5) = anrmy*rhs(2) - anrmx*rhs(6)
         b(6) = anrmy*rhs(3) - anrmx*rhs(7)
         b(7) = anrmy*rhs(4) - anrmx*rhs(8)

         ny2 = anrmy**2
         ny3 = ny2*anrmy
         ny4 = ny3*anrmy
         ny5 = ny4*anrmy

         Sx2 = (-two*b(1)*(nine - nine*ny2 + ny4) - &
              six*b(4)*anrmx*(nine - nine*ny2 + ny4) + &
              anrmy*(b(3)*anrmx*(-nine + eight*ny2) - &
              two*b(5)*(18.d0 - 26.d0*ny2 + nine*ny4) + &
              two*anrmy*(b(2)*(-one + ny2) + &
              three*b(7)*anrmy*(-one + ny2) + & 
              b(6)*anrmx*(-ten + nine*ny2))))
         
         Sy2 = -(two*(b(2) + b(6)*anrmx) + &
              (two*b(5) + six*b(7) + b(3)*anrmx)*anrmy + &
              two*(b(1) + seven*b(2) + (three*b(4) + eight*b(6))*anrmx)* &
              ny2 + two*(eight*b(5) + 21.d0*b(7) + four*b(3)*anrmx)* &
              ny3 - two*(b(1) - b(2) + &
              three*(b(4) - three*b(6))*anrmx)*ny4 + &
              six*(-three*b(5) + b(7))*ny5)

         Sxy = (b(3)*(-nine + eight*ny2)*(one + eight*ny2) + &
              two*b(5)*anrmx*(-nine + eight*ny2)*(one + nine*ny2) + &
              two*anrmy*(b(6)*(one + eight*ny2)*(-ten + nine*ny2) - &
              three*b(4)*(nine - 17.d0*ny2 + eight*ny4) + &
              anrmx*(b(1)*(-nine + eight*ny2) - &
              (b(2) + three*b(7)*anrmy)*(one + eight*ny2))))
         
         den = one / (18.d0*(-one - six*ny2 + six*ny4))
         Sx2 = Sx2 * den
         Sy2 = Sy2 * den
         Sxy = Sxy * half * den

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
       vfrac, vlo, vhi, vcent, tlo, thi, barea, alo, ahi, &
       bcent, blo, bhi, bnorm, mlo, mhi) &
       bind(c,name='amrex_eb2_build_cells')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, &
         fxlo, fxhi, fylo, fyhi, fzlo, fzhi, axlo, axhi, aylo, ayhi, azlo, azhi, &
         cxlo, cxhi, cylo, cyhi, czlo, czhi, vlo, vhi, tlo, thi, alo, ahi, blo, bhi, mlo, mhi
    integer(c_int)  , intent(inout) ::  cell( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    integer(c_int)  , intent(inout) ::    fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
    integer(c_int)  , intent(inout) ::    fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
    integer(c_int)  , intent(inout) ::    fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
    real(amrex_real), intent(inout) ::   apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(amrex_real), intent(inout) ::   apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(amrex_real), intent(inout) ::   apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(amrex_real), intent(in   ) ::   fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),3)
    real(amrex_real), intent(in   ) ::   fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),3)
    real(amrex_real), intent(in   ) ::   fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),3)
    real(amrex_real), intent(inout) :: vfrac( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3))
    real(amrex_real), intent(inout) :: vcent( tlo(1): thi(1), tlo(2): thi(2), tlo(3): thi(3),3)
    real(amrex_real), intent(inout) :: barea( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
    real(amrex_real), intent(inout) :: bcent( blo(1): bhi(1), blo(2): bhi(2), blo(3): bhi(3),3)
    real(amrex_real), intent(inout) :: bnorm( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3),3)

    integer :: i,j,k
    real(amrex_real) :: axm, axp, aym, ayp, azm, azp, aax, aay, aaz, B0, Bx, By, Bz
    real(amrex_real) :: dapx, dapy, dapz, apnorm, apnorminv, nx, ny, nz, bainv

    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             if (cell(i,j,k) .eq. regular) then
                vfrac(i,j,k) = one
                vcent(i,j,k,:) = zero
                bcent(i,j,k,:) = zero
                bnorm(i,j,k,:) = zero
                barea(i,j,k) = zero
             else if (cell(i,j,k) .eq. covered) then
                vfrac(i,j,k) = zero
                vcent(i,j,k,:) = zero
                bcent(i,j,k,:) = zero
                bnorm(i,j,k,:) = zero
                barea(i,j,k) = zero
             else
                axm = apx(i,j,k)
                axp = apx(i+1,j,k)
                aym = apy(i,j,k)
                ayp = apy(i,j+1,k)
                azm = apz(i,j,k)
                azp = apz(i,j,k+1)
                dapx = axm - axp
                dapy = aym - ayp
                dapz = azm - azp
                apnorm = sqrt(dapx**2 + dapy**2 + dapz**2)
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
                By = -ny*aay + nx*(axm*fcx(i,j,k,2)-axp*fcx(i+1,j,k,2)) &
                     &       + nz*(azm*fcz(i,j,k,2)-azp*fcz(i,j,k+1,2))
                Bz = -nz*aaz + nx*(axm*fcx(i,j,k,3)-axp*fcx(i+1,j,k,3)) &
                     &       + ny*(aym*fcy(i,j,k,3)-ayp*fcy(i,j+1,k,3))

                vfrac(i,j,k) = half*(B0 + nx*Bx + ny*By + nz*Bz)

                bainv = one / barea(i,j,k)
                bcent(i,j,k,1) = bainv * (Bx + nx*vfrac(i,j,k))
                bcent(i,j,k,2) = bainv * (By + ny*vfrac(i,j,k))
                bcent(i,j,k,3) = bainv * (Bz + nz*vfrac(i,j,k))
               
!                print *, i,j,k,vfrac(i,j,k), barea(i,j,k), bcent(i,j,k,:)

                ! xxxxx to do: vcent


                ! xxxxx clamp

                ! should we remove small cells? probably not safe to do it here
                ! because it affects face area

             end if
          end do
       end do
    end do

  end subroutine amrex_eb2_build_cells
    

end module amrex_eb2_3d_module
