module amrex_mlmg_interp_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlmg_lin_cc_interp, amrex_mlmg_lin_nd_interp

contains

  subroutine amrex_mlmg_lin_cc_interp (lo, hi, ff, fflo, ffhi, cc, cclo, cchi, ratio, nc) &
       bind(c,name='amrex_mlmg_lin_cc_interp')
    integer, dimension(1), intent(in) :: lo, hi, fflo, ffhi, cclo, cchi
    integer, intent(in) :: ratio, nc
    real(amrex_real), intent(in   ) :: cc(cclo(1):cchi(1),nc)
    real(amrex_real), intent(inout) :: ff(fflo(1):ffhi(1),nc)

    integer :: i,n, ic, ioff

    if (ratio == 2) then

       do n = 1, nc
          do i = lo(1), hi(1)
             ic = i/2
             ioff = 2*(i-ic*2)-1
             ff(i,n) = 0.75d0*cc(ic,n) + 0.25d0*cc(ic+ioff,n)
          end do
       end do
       
    else if (ratio == 4) then

       do n = 1, nc
          do i = lo(1), hi(1)
             ic = i/4
             ff(i,n) = cc(ic,n)
          end do
       end do
       
    else

       do n = 1, nc
          do i = lo(1), hi(1)
             ic = i/ratio
             ff(i,n) = cc(ic,n)
          end do
       end do

    end if

  end subroutine amrex_mlmg_lin_cc_interp


  subroutine amrex_mlmg_lin_nd_interp (clo, chi, flo, fhi, fine, fdlo, fdhi, crse, cdlo, cdhi, nc) &
       bind(c,name='amrex_mlmg_lin_nd_interp')
    integer, dimension(1) :: clo, chi, flo, fhi, fdlo, fdhi, cdlo, cdhi
    integer, intent(in) :: nc
    real(amrex_real), intent(inout) :: fine(fdlo(1):fdhi(1))
    real(amrex_real), intent(in   ) :: crse(cdlo(1):cdhi(1))
    
    integer :: i,n,ii

    do n = 1, nc
       do i = clo(1), chi(1)-1
          fine(2*i  ,n) = crse(i,n)
          fine(2*i+1,n) = 0.5d0*(crse(i,n)+crse(i+1,n))
       end do
       fine(fhi(1,n)) = crse(chi(1,n))
    end do
  end subroutine amrex_mlmg_lin_nd_interp

end module amrex_mlmg_interp_module
