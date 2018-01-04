module amrex_mlmg_interp_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlmg_lin_cc_interp, amrex_mlmg_lin_nd_interp

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


  subroutine amrex_mlmg_lin_nd_interp (clo, chi, flo, fhi, fine, fdlo, fdhi, crse, cdlo, cdhi) &
       bind(c,name='amrex_mlmg_lin_nd_interp')
    integer, dimension(2) :: clo, chi, flo, fhi, fdlo, fdhi, cdlo, cdhi
    real(amrex_real), intent(inout) :: fine(fdlo(1):fdhi(1),fdlo(2):fdhi(2))
    real(amrex_real), intent(in   ) :: crse(cdlo(1):cdhi(1),cdlo(2):cdhi(2))
    
    integer :: i,j,ii,jj

    do    j = clo(2), chi(2)
       do i = clo(1), chi(1)-1
          fine(2*i  ,2*j) = crse(i,j)
          fine(2*i+1,2*j) = 0.5d0*(crse(i,j)+crse(i+1,j))
       end do
       fine(fhi(1),2*j) = crse(chi(1),j)
       if (j < chi(2)) then
          do i = clo(1), chi(1)-1
             fine(2*i  ,2*j+1) = 0.5d0*(crse(i,j)+crse(i,j+1))
             fine(2*i+1,2*j+1) = 0.25d0*(crse(i,j)+crse(i+1,j)+crse(i,j+1)+crse(i+1,j+1))
          end do
          fine(fhi(1),2*j+1) = 0.5d0*(crse(chi(1),j)+crse(chi(1),j+1))
       end if
    end do

  end subroutine amrex_mlmg_lin_nd_interp

end module amrex_mlmg_interp_module
