module amrex_mlmg_interp_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlmg_lin_nd_interp, amrex_mlmg_eb_cc_interp

contains

  subroutine amrex_mlmg_lin_nd_interp (clo, chi, flo, fhi, fine, fdlo, fdhi, crse, cdlo, cdhi, nc) &
       bind(c,name='amrex_mlmg_lin_nd_interp')
    integer, dimension(2) :: clo, chi, flo, fhi, fdlo, fdhi, cdlo, cdhi
    integer, intent(in) :: nc
    real(amrex_real), intent(inout) :: fine(fdlo(1):fdhi(1),fdlo(2):fdhi(2),nc)
    real(amrex_real), intent(in   ) :: crse(cdlo(1):cdhi(1),cdlo(2):cdhi(2),nc)
    
    integer :: i,j,n


    ! Interpolate from coarse grid to fine grid
    do n = 1, nc
       do    j = flo(2), fhi(2)
          do i = flo(1), fhi(1)
             if (MOD(j,2) .ne. 0 .and. MOD(i,2) .ne. 0) then
                ! Fine node at center of cell
                fine(i,j,n) = 0.25d0*( &
                     crse(i/2,j/2,n) + &
                     crse(i/2 + 1, j/2, n) + &
                     crse(i/2, j/2 + 1, n) + &
                     crse(i/2 + 1, j/2 + 1, n))
             else if (MOD(j,2) .ne. 0) then
                ! Fine node along a vertical cell edge
                fine(i,j,n) = 0.5d0*( &
                     crse(i/2,j/2,n) + &
                     crse(i/2,j/2+1,n)) 
             else if (MOD(i,2) .ne. 0) then
                ! Fine node along a horizontal cell edge
                fine(i,j,n) = 0.5d0*( &
                     crse(i/2,j/2,n) + &
                     crse(i/2+1,j/2,n))
             else
                ! Fine node coincident with coarse node
                fine(i,j,n) = crse(i/2,j/2,n)
             end if
          end do
       end do
    end do

  end subroutine amrex_mlmg_lin_nd_interp


  subroutine amrex_mlmg_eb_cc_interp (lo, hi, ff, fflo, ffhi, cc, cclo, cchi, &
       flag, glo, ghi, ratio, nc) &
       bind(c,name='amrex_mlmg_eb_cc_interp')
#ifdef AMREX_USE_EB
    use amrex_ebcellflag_module, only : is_covered_cell
#endif
    integer, dimension(2), intent(in) :: lo, hi, fflo, ffhi, cclo, cchi, glo, ghi
    integer, intent(in) :: ratio, nc
    real(amrex_real), intent(in   ) :: cc(cclo(1):cchi(1),cclo(2):cchi(2),nc)
    real(amrex_real), intent(inout) :: ff(fflo(1):ffhi(1),fflo(2):ffhi(2),nc)
    integer         , intent(in   ) :: flag(glo(1):ghi(1),glo(2):ghi(2))

#ifdef AMREX_USE_EB
    
    integer :: i, j, n, ic, jc

    do n = 1, nc
       do j = lo(2), hi(2)
          jc = j/ratio
          do i = lo(1), hi(1)
             ic = i/ratio
             if (is_covered_cell(flag(i,j))) then
                ff(i,j,n) = 0.d0
             else
                ff(i,j,n) = cc(ic,jc,n)
             end if
          end do
       end do
    end do
#endif
  end subroutine amrex_mlmg_eb_cc_interp

end module amrex_mlmg_interp_module
