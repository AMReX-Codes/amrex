module amrex_mlmg_interp_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlmg_lin_cc_interp

contains

  subroutine amrex_mlmg_lin_cc_interp (lo, hi, ff, fflo, ffhi, cc, cclo, cchi, ratio) &
       bind(c,name='amrex_mlmg_lin_cc_interp')
    integer, dimension(3), intent(in) :: lo, hi, fflo, ffhi, cclo, cchi
    integer, intent(in) :: ratio
    real(amrex_real), intent(in   ) :: cc(cclo(1):cchi(1),cclo(2):cchi(2),cclo(3):cchi(3))
    real(amrex_real), intent(inout) :: ff(fflo(1):ffhi(1),fflo(2):ffhi(2),fflo(3):ffhi(3))

    integer :: i,j,k, ic, jc, kc, ioff, joff, koff

    if (ratio == 2) then

       do k = lo(3), hi(3)
          kc = k/2
          koff = 2*(k-kc*2)-1
          do j = lo(2), hi(2)
             jc = j/2
             joff = 2*(j-jc*2)-1
             do i = lo(1), hi(1)
                ic = i/2
                ioff = 2*(i-ic*2)-1
                ff(i,j,k) = 0.421875d0*cc(ic     ,jc     ,kc     ) &
                     +      0.140625d0*cc(ic+ioff,jc     ,kc     ) &
                     +      0.140625d0*cc(ic     ,jc+joff,kc     ) &
                     +      0.140625d0*cc(ic     ,jc     ,kc+koff) &
                     +      0.046875d0*cc(ic     ,jc+joff,kc+koff) &
                     +      0.046875d0*cc(ic+ioff,jc     ,kc+koff) &
                     +      0.046875d0*cc(ic+ioff,jc+joff,kc     ) &
                     +      0.015625d0*cc(ic+ioff,jc+joff,kc+koff)
             end do
          end do
       end do
       
    else

       do k = lo(3), hi(3)
          kc = k/2
          do j = lo(2), hi(2)
             jc = j/2
             do i = lo(1), hi(1)
                ic = i/2
                ff(i,j,k) = cc(ic,jc,kc)
             end do
          end do
       end do
       
    end if

  end subroutine amrex_mlmg_lin_cc_interp

end module amrex_mlmg_interp_module
