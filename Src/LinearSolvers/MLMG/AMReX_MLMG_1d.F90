module amrex_mlmg_interp_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlmg_lin_nd_interp, amrex_mlmg_eb_cc_interp

contains

  subroutine amrex_mlmg_lin_nd_interp (clo, chi, flo, fhi, fine, fdlo, fdhi, crse, cdlo, cdhi, nc) &
       bind(c,name='amrex_mlmg_lin_nd_interp')
    integer, dimension(1) :: clo, chi, flo, fhi, fdlo, fdhi, cdlo, cdhi
    integer, intent(in) :: nc
    real(amrex_real), intent(inout) :: fine(fdlo(1):fdhi(1),nc)
    real(amrex_real), intent(in   ) :: crse(cdlo(1):cdhi(1),nc)
    
    integer :: i,n,ii

    do n = 1, nc
       do i = clo(1), chi(1)-1
          fine(2*i  ,n) = crse(i,n)
          fine(2*i+1,n) = 0.5d0*(crse(i,n)+crse(i+1,n))
       end do
       fine(fhi(1),n) = crse(chi(1),n)
    end do
  end subroutine amrex_mlmg_lin_nd_interp


  subroutine amrex_mlmg_eb_cc_interp (lo, hi, ff, fflo, ffhi, cc, cclo, cchi, &
       flag, glo, ghi, ratio, nc) &
       bind(c,name='amrex_mlmg_eb_cc_interp')
    integer, dimension(1), intent(in) :: lo, hi, fflo, ffhi, cclo, cchi, glo, ghi
    integer, intent(in) :: ratio, nc
    real(amrex_real), intent(in   ) :: cc(cclo(1):cchi(1),nc)
    real(amrex_real), intent(inout) :: ff(fflo(1):ffhi(1),nc)
    integer, intent(in) :: flag(glo(1):ghi(1))

    call amrex_error("1d eb not supported")
  end subroutine amrex_mlmg_eb_cc_interp

end module amrex_mlmg_interp_module
