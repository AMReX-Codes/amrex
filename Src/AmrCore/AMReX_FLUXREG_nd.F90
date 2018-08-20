module amrex_fluxreg_nd_module

  use amrex_fort_module
  implicit none

  integer, parameter, private :: crse_cell = 0
  integer, parameter, private :: fine_ontop = 1

contains

  subroutine amrex_froverwrite_cfb (lo, hi, dst, dlo, dhi, src, slo, shi, &
       msk, mlo, mhi, nc, dir, scale) bind(c,name='amrex_froverwrite_cfb')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, slo, shi, mlo, mhi
    integer, intent(in) :: nc, dir
    real(amrex_real), intent(in) :: scale
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),nc)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nc)
    integer,          intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k,n

    if (dir .eq. 0) then
       do n = 1, nc
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if ( (msk(i-1,j,k).eq.crse_cell  .and. msk(i,j,k).eq.fine_ontop) .or. &
                        (msk(i-1,j,k).eq.fine_ontop .and. msk(i,j,k).eq.crse_cell ) ) then
                      dst(i,j,k,n) = scale*src(i,j,k,n)
                   end if
                end do
             end do
          end do
       end do
    else if (dir .eq. 1) then
       do n = 1, nc
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if ( (msk(i,j-1,k).eq.crse_cell  .and. msk(i,j,k).eq.fine_ontop) .or. &
                        (msk(i,j-1,k).eq.fine_ontop .and. msk(i,j,k).eq.crse_cell ) ) then
                      dst(i,j,k,n) = scale*src(i,j,k,n)
                   end if
                end do
             end do
          end do
       end do
    else
       do n = 1, nc
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   if ( (msk(i,j,k-1).eq.crse_cell  .and. msk(i,j,k).eq.fine_ontop) .or. &
                        (msk(i,j,k-1).eq.fine_ontop .and. msk(i,j,k).eq.crse_cell ) ) then
                      dst(i,j,k,n) = scale*src(i,j,k,n)
                   end if
                end do
             end do
          end do
       end do
    end if
  end subroutine amrex_froverwrite_cfb

end module amrex_fluxreg_nd_module
