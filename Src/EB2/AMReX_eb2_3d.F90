
module amrex_eb2_3d_module

  use amrex_fort_module
  implicit none

  integer, parameter :: regular = 0
  integer, parameter :: covered = 1
  integer, parameter :: irregular = 2
  integer, parameter :: unknown = 3

  private
  public :: amrex_eb2_gfab_build_types

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

end module amrex_eb2_3d_module
