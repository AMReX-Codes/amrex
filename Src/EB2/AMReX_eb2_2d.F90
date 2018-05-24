
module amrex_eb2_2d_moudle

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
       fx, fxlo, fxhi, fy, fylo, fyhi) bind(c,name='amrex_eb2_gfab_build_types')
    integer, dimension(2), intent(in) :: lo, hi, slo, shi, clo, chi, fxlo, fxhi, fylo, fyhi
    real(amrex_real), intent(in   ) ::    s( slo(1): shi(1), slo(2): shi(2))
    integer(c_int)  , intent(inout) :: cell( clo(1): chi(1), clo(2): chi(2))
    integer(c_int)  , intent(inout) ::   fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    integer(c_int)  , intent(inout) ::   fy(fylo(1):fyhi(1),fylo(2):fyhi(2))

    integer :: i,j

    do    j = lo(2)-2, hi(2)+2
       do i = lo(1)-2, hi(1)+2
          if (       s(i,j  ).ge.0.d0 .and. s(i+1,j  ).ge.0.d0 &
               .and. s(i,j+1).ge.0.d0 .and. s(i+1,j+1).ge.0.d0) then
             cell(i,j) = covered
          else if (  s(i,j  ).lt.0.d0 .and. s(i+1,j  ).lt.0.d0 &
               .and. s(i,j+1).lt.0.d0 .and. s(i+1,j+1).lt.0.d0) then
             cell(i,j) = regular
          else
             cell(i,j) = irregular
          end if
       end do
    end do

    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)-1, hi(1)+2
          if (s(i,j).ge.0.d0 .and. s(i,j+1).ge.0.d0) then
             fx(i,j) = covered
          else if (s(i,j).lt.0.d0 .and. s(i,j+1).lt.0.d0) then
             fx(i,j) = regular
          else
             fx(i,j) = irregular
          end if
       end do
    end do

    do    j = lo(2)-1, hi(2)+2
       do i = lo(1)-1, hi(1)+1
          if (s(i,j).ge.0.d0 .and. s(i+1,j).ge.0.d0) then
             fy(i,j) = covered
          else if (s(i,j).lt.0.d0 .and. s(i+1,j).lt.0.d0) then
             fy(i,j) = regular
          else
             fy(i,j) = irregular
          end if
       end do
    end do

  end subroutine amrex_eb2_gfab_build_types

end module amrex_eb2_2d_moudle
