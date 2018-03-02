module amrex_mlmg_interp_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlmg_lin_cc_interp, amrex_mlmg_lin_nd_interp

contains

  subroutine amrex_mlmg_lin_cc_interp (lo, hi, &
       ff, fflo, ffhi, &
       cc, cclo, cchi, &
       ratio, nc) &
       bind(c,name='amrex_mlmg_lin_cc_interp')
    integer, dimension(2), intent(in) :: lo, hi, fflo, ffhi, cclo, cchi
    integer, intent(in) :: ratio, nc
    real(amrex_real), intent(in   ) :: cc(cclo(1):cchi(1),cclo(2):cchi(2),nc)
    real(amrex_real), intent(inout) :: ff(fflo(1):ffhi(1),fflo(2):ffhi(2),nc)

    integer :: i,j,n, ic, jc, ioff, joff

    if (ratio == 2) then
       do n = 1, nc
          do j = lo(2), hi(2)
             jc = j/2
             joff = 2*(j-jc*2)-1
             do i = lo(1), hi(1)
                ic = i/2
                ioff = 2*(i-ic*2)-1
                ff(i,j,n) = 0.5625d0*cc(ic     ,jc     ,n) &
                       +    0.1875d0*cc(ic+ioff,jc     ,n) &
                       +    0.1875d0*cc(ic     ,jc+joff,n) &
                       +    0.0625d0*cc(ic+ioff,jc+joff,n)
             end do
          end do
       end do
       
    else if (ratio == 4) then

       do n = 1, 1
          do j = lo(2), hi(2)
             jc = j/4
             do i = lo(1), hi(1)
                ic = i/4
                ff(i,j,n) = cc(ic,jc,n)
             end do
          end do
       end do
       
    else

       call amrex_abort("amrex_mlmg_lin_cc_interp: only ratio 2 and 4 are supported")

    end if

  end subroutine amrex_mlmg_lin_cc_interp


  subroutine amrex_mlmg_lin_nd_interp (clo, chi, flo, fhi, fine, fdlo, fdhi, crse, cdlo, cdhi, nc) &
       bind(c,name='amrex_mlmg_lin_nd_interp')
    integer, dimension(2) :: clo, chi, flo, fhi, fdlo, fdhi, cdlo, cdhi
    integer, intent(in) :: nc
    real(amrex_real), intent(inout) :: fine(fdlo(1):fdhi(1),fdlo(2):fdhi(2),nc)
    real(amrex_real), intent(in   ) :: crse(cdlo(1):cdhi(1),cdlo(2):cdhi(2),nc)
    
    integer :: i,j,n,ii,jj

    do n = 1, nc
       do    j = clo(2), chi(2)
          do i = clo(1), chi(1)-1
             fine(2*i  ,2*j,n) = crse(i,j,n)
             fine(2*i+1,2*j,n) = 0.5d0*(crse(i,j,n)+crse(i+1,j,n))
          end do
          fine(fhi(1),2*j,n) = crse(chi(1),j,n)
          if (j < chi(2)) then
             do i = clo(1), chi(1)-1
                fine(2*i  ,2*j+1,n) = 0.5d0*(crse(i,j,n)+crse(i,j+1,n))
                fine(2*i+1,2*j+1,n) = 0.25d0*(crse(i,j,n)+crse(i+1,j,n)+crse(i,j+1,n)+crse(i+1,j+1,n))
             end do
             fine(fhi(1),2*j+1,n) = 0.5d0*(crse(chi(1),j,n)+crse(chi(1),j+1,n))
          end if
       end do
    end do
    
  end subroutine amrex_mlmg_lin_nd_interp

end module amrex_mlmg_interp_module
