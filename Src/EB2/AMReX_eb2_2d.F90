
module amrex_eb2_2d_moudle

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module, only : zero, one, half
  implicit none

  integer, parameter :: regular = 0
  integer, parameter :: covered = 1
  integer, parameter :: irregular = 2
  integer, parameter :: unknown = 3

  real(amrex_real), private, parameter :: small = 1.d-14

  private
  public :: amrex_eb2_gfab_build_types, amrex_eb2_build_faces

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

end module amrex_eb2_2d_moudle
