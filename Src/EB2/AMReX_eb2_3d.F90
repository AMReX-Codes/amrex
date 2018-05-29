
module amrex_eb2_3d_module

  use amrex_error_module
  use amrex_fort_module
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
             if (       s(i,j  ,k  ).ge.0.d0 .and. s(i+1,j  ,k  ).ge.0.d0 &
                  .and. s(i,j+1,k  ).ge.0.d0 .and. s(i+1,j+1,k  ).ge.0.d0 &
                  .and. s(i,j  ,k+1).ge.0.d0 .and. s(i+1,j  ,k+1).ge.0.d0 &
                  .and. s(i,j+1,k+1).ge.0.d0 .and. s(i+1,j+1,k+1).ge.0.d0 ) then
                cell(i,j,k) = covered
             else if (  s(i,j  ,k  ).lt.0.d0 .and. s(i+1,j  ,k  ).lt.0.d0 &
                  .and. s(i,j+1,k  ).lt.0.d0 .and. s(i+1,j+1,k  ).lt.0.d0 &
                  .and. s(i,j  ,k+1).lt.0.d0 .and. s(i+1,j  ,k+1).lt.0.d0 &
                  .and. s(i,j+1,k+1).lt.0.d0 .and. s(i+1,j+1,k+1).lt.0.d0 ) then
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
             if (       s(i,j,k  ).ge.0.d0 .and. s(i,j+1,k  ).ge.0.d0 &
                  .and. s(i,j,k+1).ge.0.d0 .and. s(i,j+1,k+1).ge.0.d0 ) then
                fx(i,j,k) = covered
             else if (  s(i,j,k  ).lt.0.d0 .and. s(i,j+1,k  ).lt.0.d0 &
                  .and. s(i,j,k+1).lt.0.d0 .and. s(i,j+1,k+1).lt.0.d0 ) then
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
             if (       s(i,j,k  ).ge.0.d0 .and. s(i+1,j,k  ).ge.0.d0 &
                  .and. s(i,j,k+1).ge.0.d0 .and. s(i+1,j,k+1).ge.0.d0 ) then
                fy(i,j,k) = covered
             else if (  s(i,j,k  ).lt.0.d0 .and. s(i+1,j,k  ).lt.0.d0 &
                  .and. s(i,j,k+1).lt.0.d0 .and. s(i+1,j,k+1).lt.0.d0 ) then
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
             if (       s(i,j  ,k).ge.0.d0 .and. s(i+1,j  ,k).ge.0.d0 &
                  .and. s(i,j+1,k).ge.0.d0 .and. s(i+1,j+1,k).ge.0.d0 ) then
                fz(i,j,k) = covered
             else if (  s(i,j  ,k).lt.0.d0 .and. s(i+1,j  ,k).lt.0.d0 &
                  .and. s(i,j+1,k).lt.0.d0 .and. s(i+1,j+1,k).lt.0.d0 ) then
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
             if (s(i,j,k).ge.0.d0 .and. s(i+1,j,k).ge.0.d0) then
                ex(i,j,k) = covered
             else if (s(i,j,k).lt.0.d0 .and. s(i+1,j,k).lt.0.d0) then
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
             if (s(i,j,k).ge.0.d0 .and. s(i,j+1,k).ge.0.d0) then
                ey(i,j,k) = covered
             else if (s(i,j,k).lt.0.d0 .and. s(i,j+1,k).lt.0.d0) then
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
             if (s(i,j,k).ge.0.d0 .and. s(i,j,k+1).ge.0.d0) then
                ez(i,j,k) = covered
             else if (s(i,j,k).lt.0.d0 .and. s(i,j,k+1).lt.0.d0) then
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
       dx, dxinv, problo) &
       bind(c,name='amrex_eb2_build_faces')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, slo, shi, &
         fxlo, fxhi, fylo, fyhi, fzlo, fzhi, exlo, exhi, eylo, eyhi, ezlo, ezhi, &
         ixlo, ixhi, iylo, iyhi, izlo, izhi, axlo, axhi, aylo, ayhi, azlo, azhi, &
         cxlo, cxhi, cylo, cyhi, czlo, czhi
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

    integer :: i,j,k, ncuts
    real(amrex_real) :: cut
    real(amrex_real) :: lxm, lxp, lym, lyp, lzm, lzp, bcx, bcy, bcz

    ! x-face
    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+2
             if (fx(i,j,k) .eq. regular) then
                apx(i,j,k) = 1.d0
                fcx(i,j,k,:) = 0.d0
             else if (fx(i,j,k) .eq. covered) then
                apx(i,j,k) = 0.d0
                fcx(i,j,k,:) = 0.d0
             else
                ncuts = 0
                bcy = 0.d0
                bcz = 0.d0

                if (ey(i,j,k) .eq. regular) then
                   lym = 1.d0
                else if (ey(i,j,k) .eq. covered) then
                   lym = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (intery(i,j,k)-(problo(2)+j*dx(2)))*dxinv(2)
                   bcy = bcy + cut
                   if (levset(i,j,k) .lt. 0.d0) then
                      lym = cut
                   else
                      lym = 1.d0-cut
                   end if
                   lym = min(max(0.d0,lym),1.d0)
                end if

                if (ey(i,j,k+1) .eq. regular) then
                   lyp = 1.d0
                else if (ey(i,j,k+1) .eq. covered) then
                   lyp = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (intery(i,j,k+1)-(problo(2)+j*dx(2)))*dxinv(2)
                   bcy = bcy + cut
                   bcz = bcz + 1.d0
                   if (levset(i,j,k+1) .lt. 0.d0) then
                      lyp = cut
                   else
                      lyp = 1.d0-cut
                   end if
                   lyp = min(max(0.d0,lyp),1.d0)
                end if

                if (ez(i,j,k) .eq. regular) then
                   lzm = 1.d0
                else if (ez(i,j,k) .eq. covered) then
                   lzm = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (interz(i,j,k)-(problo(3)+k*dx(3)))*dxinv(3)
                   bcz = bcz + cut
                   if (levset(i,j,k) .lt. 0.d0) then
                      lzm = cut
                   else
                      lzm = 1.d0-cut
                   end if
                   lzm = min(max(0.d0,lzm),1.d0)
                end if

                if (ez(i,j+1,k) .eq. regular) then
                   lzp = 1.d0
                else if (ez(i,j+1,k) .eq. covered) then
                   lzp = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (interz(i,j+1,k)-(problo(3)+k*dx(3)))*dxinv(3)
                   bcy = bcy + 1.d0
                   bcz = bcz + cut
                   if (levset(i,j+1,k) .lt. 0.d0) then
                      lzp = cut
                   else
                      lzp = 1.d0-cut
                   end if
                end if

                if (ncuts .gt. 2) then
                   call amrex_error("amrex_eb2_build_faces: more than 2 cuts not suported")
                else if (ncuts .lt. 2) then
                   print *, "ncuts = ", i,j,k,ncuts
                   call amrex_error("amrex_eb2_build_faces: irregular face with less than 2 cuts???")
                end if

                if (lym.eq.lyp .and. lzm.eq.lzp) then
                   apx(i,j,k) = 1.d0
                   fcx(i,j,k,:) = 0.d0
                else
                   fcx(i,j,k,1) = 0.d0
                   bcy = 0.5d0*bcy - 0.5d0
                   bcz = 0.5d0*bcz - 0.5d0
                   call cut_face_2d(apx(i,j,k),fcx(i,j,k,2),fcx(i,j,k,3),lzm,lzp,lym,lyp,bcy,bcz)
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
                apy(i,j,k) = 1.d0
                fcy(i,j,k,:) = 0.d0
             else if (fy(i,j,k) .eq. covered) then
                apy(i,j,k) = 0.d0
                fcy(i,j,k,:) = 0.d0
             else
                ncuts = 0
                bcx = 0.d0
                bcz = 0.d0

                if (ex(i,j,k) .eq. regular) then
                   lxm = 1.d0
                else if (ex(i,j,k) .eq. covered) then
                   lxm = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (interx(i,j,k)-(problo(1)+i*dx(1)))*dxinv(1)
                   bcx = bcx + cut
                   if (levset(i,j,k) .lt. 0.d0) then
                      lxm = cut
                   else
                      lxm = 1.d0-cut
                   end if
                   lxm = min(max(0.d0,lxm),1.d0)
                end if

                if (ex(i,j,k+1) .eq. regular) then
                   lxp = 1.d0
                else if (ex(i,j,k+1) .eq. covered) then
                   lxp = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (interx(i,j,k+1)-(problo(1)+i*dx(1)))*dxinv(1)
                   bcx = bcx + cut
                   bcz = bcz + 1.d0
                   if (levset(i,j,k+1) .lt. 0.d0) then
                      lxp = cut
                   else
                      lxp = 1.d0-cut
                   end if
                   lxp = min(max(0.d0,lxp),1.d0)
                end if

                if (ez(i,j,k) .eq. regular) then
                   lzm = 1.d0
                else if (ez(i,j,k) .eq. covered) then
                   lzm = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (interz(i,j,k)-(problo(3)+k*dx(3)))*dxinv(3)
                   bcz = bcz + cut
                   if (levset(i,j,k) .lt. 0.d0) then
                      lzm = cut
                   else
                      lzm = 1.d0-cut
                   end if
                end if

                if (ez(i+1,j,k) .eq. regular) then
                   lzp = 1.d0
                else if (ez(i+1,j,k) .eq. covered) then
                   lzp = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (interz(i+1,j,k)-(problo(3)+k*dx(3)))*dxinv(3)
                   bcx = bcx + 1.d0
                   bcz = bcz + cut
                   if (levset(i+1,j,k) .lt. 0.d0) then
                      lzp = cut
                   else
                      lzp = 1.d0-cut
                   end if
                end if

                if (ncuts .gt. 2) then
                   call amrex_error("amrex_eb2_build_faces: more than 2 cuts not supported")
                else if (ncuts .lt. 2) then
                   print *, "ncuts = ", i,j,k,ncuts
                   call amrex_error("amrex_eb2_build_faces: irregular face with less than 2 cuts???")
                end if

                if (lxm.eq.lxp .and. lzm.eq.lzp) then
                   apy(i,j,k) = 1.d0
                   fcy(i,j,k,:) = 0.d0
                else
                   fcy(i,j,k,2) = 0.d0
                   bcx = 0.5d0*bcx - 0.5d0
                   bcz = 0.5d0*bcz - 0.5d0
                   call cut_face_2d(apy(i,j,k),fcy(i,j,k,1),fcy(i,j,k,3),lzm,lzp,lxm,lxp,bcx,bcz)
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
                apz(i,j,k) = 1.d0
                fcz(i,j,k,:) = 0.d0
             else if (fz(i,j,k) .eq. covered) then
                apz(i,j,k) = 0.d0
                fcz(i,j,k,:) = 0.d0
             else
                ncuts = 0
                bcx = 0.d0
                bcy = 0.d0

                if (ex(i,j,k) .eq. regular) then
                   lxm = 1.d0
                else if (ex(i,j,k) .eq. covered) then
                   lxm = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (interx(i,j,k)-(problo(1)+i*dx(1)))*dxinv(1)
                   bcx = bcx + cut
                   if (levset(i,j,k) .lt. 0.d0) then
                      lxm = cut
                   else
                      lxm = 1.d0-cut
                   end if
                   lxm = min(max(0.d0,lxm),1.d0)
                end if

                if (ex(i,j+1,k) .eq. regular) then
                   lxp = 1.d0
                else if (ex(i,j+1,k) .eq. covered) then
                   lxp = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (interx(i,j+1,k)-(problo(1)+i*dx(1)))*dxinv(1)
                   bcx = bcx + cut
                   bcy = bcy + 1.d0
                   if (levset(i,j+1,k) .lt. 0.d0) then
                      lxp = cut
                   else
                      lxp = 1.d0-cut
                   end if
                   lxp = min(max(0.d0,lxp),1.d0)
                end if

                if (ey(i,j,k) .eq. regular) then
                   lym = 1.d0
                else if (ey(i,j,k) .eq. covered) then
                   lym = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (intery(i,j,k)-(problo(2)+j*dx(2)))*dxinv(2)
                   bcy = bcy + cut
                   if (levset(i,j,k) .lt. 0.d0) then
                      lym = cut
                   else
                      lym = 1.d0-cut
                   end if
                end if

                if (ey(i+1,j,k) .eq. regular) then
                   lyp = 1.d0
                else if (ey(i+1,j,k) .eq. covered) then
                   lyp = 0.d0
                else
                   ncuts = ncuts+1
                   cut = (intery(i+1,j,k)-(problo(2)+j*dx(2)))*dxinv(2)
                   bcx = bcx + 1.d0
                   bcy = bcy + cut
                   if (levset(i+1,j,k) .lt. 0.d0) then
                      lyp = cut
                   else
                      lyp = 1.d0-cut
                   end if
                end if

                if (ncuts .gt. 2) then
                   call amrex_error("amrex_eb2_build_faces: more than 2 cuts not supported")
                else if (ncuts .lt. 2) then
                   print *, "ncuts = ", i,j,k,ncuts
                   call amrex_error("amrex_eb2_build_faces: irregular face with less than 2 cuts???")
                end if

                if (lxm.eq.lxp .and. lym.eq.lyp) then
                   apz(i,j,k) = 1.d0
                   fcz(i,j,k,:) = 0.d0
                else
                   fcz(i,j,k,3) = 0.d0
                   bcx = 0.5d0*bcx - 0.5d0
                   bcy = 0.5d0*bcy - 0.5d0
                   call cut_face_2d(apz(i,j,k),fcz(i,j,k,1),fcz(i,j,k,2),lym,lyp,lxm,lxp,bcx,bcy)
                end if
             end if
          end do
       end do
    end do

    do       k = lo(3)-1, hi(3)+1
       do    j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1
             if (cell(i,j,k) .eq. irregular) then
                if  (apx(i,j,k).eq.0.d0 .and. apx(i+1,j,k).eq.0.d0 .and. &
                     apy(i,j,k).eq.0.d0 .and. apy(i,j+1,k).eq.0.d0 .and. &
                     apz(i,j,k).eq.0.d0 .and. apz(i,j,k+1).eq.0.d0) then
                   cell(i,j,k) = covered
                   fx(i,j,k) = covered
                   fx(i+1,j,k) = covered
                   fy(i,j,k) = covered
                   fy(i,j+1,k) = covered
                   fz(i,j,k) = covered
                   fz(i,j,k+1) = covered
                else if (apx(i,j,k).eq.1.d0 .and. apx(i+1,j,k).eq.1.d0 .and. &
                     &   apy(i,j,k).eq.1.d0 .and. apy(i,j+1,k).eq.1.d0 .and. &
                     &   apz(i,j,k).eq.1.d0 .and. apz(i,j,k+1).eq.1.d0) then
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

    subroutine cut_face_2d (areafrac, centx, centy, axm, axp, aym, ayp, bcx, bcy)
      real(amrex_real), intent(out) :: areafrac, centx, centy
      real(amrex_real), intent(in) :: axm, axp, aym, ayp, bcx, bcy

      real(amrex_real) :: apnorm, apnorminv, anrmx, anrmy, x_ym, x_yp, y_xm, y_xp, aa, af1, af2

      apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2)
      apnorminv = 1.d0/apnorm
      anrmx = (axm-axp) * apnorminv  ! pointing to the wall
      anrmy = (aym-ayp) * apnorminv

      if (anrmx .eq. 0.d0) then
         areafrac = axm
         centx = 0.d0
         centy = (0.125d0*(ayp-aym) + anrmy*0.5d0*bcy**2)/areafrac
      else if (anrmy .eq. 0.d0) then
         areafrac = aym
         centx = (0.125d0*(axp-axm) + anrmx*0.5d0*bcx**2)/areafrac
         centy = 0.d0
      else
         if (anrmx .gt. 0.d0) then
            x_ym = -0.5d0 + aym
            x_yp = -0.5d0 + ayp
         else
            x_ym = 0.5d0 - aym
            x_yp = 0.5d0 - ayp
         end if
         aa = abs(anrmx)/anrmy
         af1 = 0.5d0*(axm+axp) + aa*0.5d0*(x_ym**2-x_yp**2)
         centx = 0.125d0*(axp-axm) + aa*(1.d0/6.d0)*(x_ym**3-x_yp**3)

         if (anrmy .gt. 0.d0) then
            y_xm = -0.5d0 + axm
            y_xp = -0.5d0 + axp
         else
            y_xm = 0.5d0 - axm
            y_xp = 0.5d0 - axp
         end if
         aa = abs(anrmy)/anrmx
         af2 = 0.5d0*(aym+ayp) + aa*0.5d0*(y_xm**2-y_xp**2)
         centy = 0.125d0*(ayp-aym) + aa*(1.d0/6.d0)*(y_xm**3-y_xp**3)

         areafrac = 0.5d0*(af1+af2)
         if (areafrac .gt. 1.d0-small) then
            areafrac = 1.d0
            centx = 0.d0
            centy = 0.d0
         else if (areafrac .lt. small) then
            areafrac = 0.d0
            centx = 0.d0
            centy = 0.d0
         else
            centx = centx*(1.d0/areafrac)
            centy = centy*(1.d0/areafrac)
            centx = min(max(centx,-0.5d0),0.5d0)
            centy = min(max(centy,-0.5d0),0.5d0)
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
                vfrac(i,j,k) = 1.d0
                vcent(i,j,k,:) = 0.d0
                bcent(i,j,k,:) = 0.d0
                bnorm(i,j,k,:) = 0.d0
                barea(i,j,k) = 0.d0
             else if (cell(i,j,k) .eq. covered) then
                vfrac(i,j,k) = 0.d0
                vcent(i,j,k,:) = 0.d0
                bcent(i,j,k,:) = 0.d0
                bnorm(i,j,k,:) = 0.d0
                barea(i,j,k) = 0.d0
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
                apnorminv = 1.d0/apnorm
                nx = dapx * apnorminv
                ny = dapy * apnorminv
                nz = dapz * apnorminv
                bnorm(i,j,k,1) = nx
                bnorm(i,j,k,2) = ny
                bnorm(i,j,k,3) = nz
                barea(i,j,k) = nx*dapx + ny*dapy + nz*dapz

                aax = 0.5d0*(axm+axp)
                aay = 0.5d0*(aym+ayp)
                aaz = 0.5d0*(azm+azp)
                B0 = aax + aay + aaz
                Bx = -nx*aax + ny*(aym*fcy(i,j,k,1)-ayp*fcy(i,j+1,k,1)) &
                     &       + nz*(azm*fcz(i,j,k,1)-azp*fcz(i,j,k+1,1))
                By = -ny*aay + nx*(axm*fcx(i,j,k,2)-axp*fcx(i+1,j,k,2)) &
                     &       + nz*(azm*fcz(i,j,k,2)-azp*fcz(i,j,k+1,2))
                Bz = -nz*aaz + nx*(axm*fcx(i,j,k,3)-axp*fcx(i+1,j,k,3)) &
                     &       + ny*(aym*fcy(i,j,k,3)-ayp*fcy(i,j+1,k,3))

                vfrac(i,j,k) = 0.5d0*(B0 + nx*Bx + ny*By + nz*Bz)

                bainv = 1.d0 / barea(i,j,k)
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
