module amrex_mlmg_interp_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlmg_lin_cc_interp

contains

  subroutine amrex_mlmg_lin_cc_interp (lo, hi, ff, fflo, ffhi, cc, cclo, cchi, ratio) &
       bind(c,name='amrex_mlmg_lin_cc_interp')
    integer, dimension(1), intent(in) :: lo, hi, fflo, ffhi, cclo, cchi
    integer, intent(in) :: ratio
    real(amrex_real), intent(in   ) :: cc(cclo(1):cchi(1))
    real(amrex_real), intent(inout) :: ff(fflo(1):ffhi(1))

    integer :: i, ic, ioff

    if (ratio == 2) then

       do i = lo(1), hi(1)
          ic = i/2
          ioff = 2*(i-ic*2)-1
          ff(i) = 0.75d0*cc(ic) + 0.25d0*cc(ic+ioff)
       end do
       
    else if (ratio == 4) then

       do i = lo(1), hi(1)
          ic = i/4
          ff(i) = cc(ic)
       end do
       
    else

       do i = lo(1), hi(1)
          ic = i/ratio
          ff(i) = cc(ic)
       end do

    end if

  end subroutine amrex_mlmg_lin_cc_interp

end module amrex_mlmg_interp_module
