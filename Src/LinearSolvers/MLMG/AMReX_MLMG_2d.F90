module amrex_mlmg_interp_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlmg_lin_cc_interp

contains

  subroutine amrex_mlmg_lin_cc_interp (lo, hi, ff, fflo, ffhi, cc, cclo, cchi, ratio) &
       bind(c,name='amrex_mlmg_lin_cc_interp')
    integer, dimension(2), intent(in) :: lo, hi, fflo, ffhi, cclo, cchi
    integer, intent(in) :: ratio
    real(amrex_real), intent(in   ) :: cc(cclo(1):cchi(1),cclo(2):cchi(2))
    real(amrex_real), intent(inout) :: ff(fflo(1):ffhi(1),fflo(2):ffhi(2))

    integer :: i,j, ic, jc, ioff, joff

    if (ratio == 2) then

       do j = lo(2), hi(2)
          jc = j/2
          joff = 2*(j-jc*2)-1
          do i = lo(1), hi(1)
             ic = i/2
             ioff = 2*(i-ic*2)-1
             ff(i,j) = 0.5625d0*cc(ic     ,jc     ) &
                  +    0.1875d0*cc(ic+ioff,jc     ) &
                  +    0.1875d0*cc(ic     ,jc+joff) &
                  +    0.0625d0*cc(ic+ioff,jc+joff)
          end do
       end do
       
    else if (ratio == 4) then

       do j = lo(2), hi(2)
          jc = j/4
          do i = lo(1), hi(1)
             ic = i/4
             ff(i,j) = cc(ic,jc)
          end do
       end do
       
    else

       call amrex_abort("amrex_mlmg_lin_cc_interp: only ratio 2 and 4 are supported")

    end if

  end subroutine amrex_mlmg_lin_cc_interp

end module amrex_mlmg_interp_module
