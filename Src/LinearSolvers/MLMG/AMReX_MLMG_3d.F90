module amrex_mlmg_interp_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlmg_lin_cc_interp, amrex_mlmg_lin_nd_interp

contains

  subroutine amrex_mlmg_lin_cc_interp (lo, hi, ff, fflo, ffhi, cc, cclo, cchi, ratio, nc) &
       bind(c,name='amrex_mlmg_lin_cc_interp')
    integer, dimension(3), intent(in) :: lo, hi, fflo, ffhi, cclo, cchi
    integer, intent(in) :: ratio, nc
    real(amrex_real), intent(in   ) :: cc(cclo(1):cchi(1),cclo(2):cchi(2),cclo(3):cchi(3),nc)
    real(amrex_real), intent(inout) :: ff(fflo(1):ffhi(1),fflo(2):ffhi(2),fflo(3):ffhi(3),nc)

    integer :: i,j,k,n, ic, jc, kc, ioff, joff, koff

    if (ratio == 2) then

       do n = 1, nc
          do k = lo(3), hi(3)
             kc = k/2
             koff = 2*(k-kc*2)-1
             do j = lo(2), hi(2)
                jc = j/2
                joff = 2*(j-jc*2)-1
                do i = lo(1), hi(1)
                   ic = i/2
                   ioff = 2*(i-ic*2)-1
                   ff(i,j,k,n) = 0.421875d0*cc(ic     ,jc     ,kc     ,n) &
                          +      0.140625d0*cc(ic+ioff,jc     ,kc     ,n) &
                          +      0.140625d0*cc(ic     ,jc+joff,kc     ,n) &
                          +      0.140625d0*cc(ic     ,jc     ,kc+koff,n) &
                          +      0.046875d0*cc(ic     ,jc+joff,kc+koff,n) &
                          +      0.046875d0*cc(ic+ioff,jc     ,kc+koff,n) &
                          +      0.046875d0*cc(ic+ioff,jc+joff,kc     ,n) &
                          +      0.015625d0*cc(ic+ioff,jc+joff,kc+koff,n)
                end do
             end do
          end do
       end do

    else if (ratio == 4) then

       do n = 1, nc
          do k = lo(3), hi(3)
             kc = k/4
             do j = lo(2), hi(2)
                jc = j/4
                do i = lo(1), hi(1)
                   ic = i/4
                   ff(i,j,k,n) = cc(ic,jc,kc,n)
                end do
             end do
          end do
       end do

    else

       call amrex_abort("amrex_mlmg_lin_cc_interp: only ratio 2 and 4 are supported")

    end if

  end subroutine amrex_mlmg_lin_cc_interp


  subroutine amrex_mlmg_lin_nd_interp (clo, chi, flo, fhi, fine, fdlo, fdhi, crse, cdlo, cdhi, nc) &
       bind(c,name='amrex_mlmg_lin_nd_interp')
    integer, dimension(3) :: clo, chi, flo, fhi, fdlo, fdhi, cdlo, cdhi
    integer, intent(in) :: nc
    real(amrex_real), intent(inout) :: fine(fdlo(1):fdhi(1),fdlo(2):fdhi(2),fdlo(3):fdhi(3),nc)
    real(amrex_real), intent(in   ) :: crse(cdlo(1):cdhi(1),cdlo(2):cdhi(2),cdlo(3):cdhi(3),nc)
    
    integer :: i,j,k,n,ii,jj,kk

    do n = 1, nc
       do k = clo(3), chi(3)-1
          kk = k*2
          do j = clo(2), chi(2)-1
             jj = j*2
             do i = clo(1), chi(1)-1
                ii = i*2
                fine(ii  ,jj  ,kk  ,n) = crse(i,j,k,n)
                fine(ii+1,jj  ,kk  ,n) = 0.5d0  *(crse(i,j  ,k  ,n)+crse(i+1,j  ,k  ,n))
                fine(ii  ,jj+1,kk  ,n) = 0.5d0  *(crse(i,j  ,k  ,n)+crse(i  ,j+1,k  ,n))
                fine(ii+1,jj+1,kk  ,n) = 0.25d0 *(crse(i,j  ,k  ,n)+crse(i+1,j  ,k  ,n) &
                     &                         +crse(i,j+1,k  ,n)+crse(i+1,j+1,k  ,n))
                fine(ii  ,jj  ,kk+1,n) = 0.5d0  *(crse(i,j  ,k  ,n)+crse(i  ,j  ,k+1,n))
                fine(ii+1,jj  ,kk+1,n) = 0.25d0 *(crse(i,j  ,k  ,n)+crse(i+1,j  ,k  ,n) &
                     &                         +crse(i,j  ,k+1,n)+crse(i+1,j  ,k+1,n))
                fine(ii  ,jj+1,kk+1,n) = 0.25d0 *(crse(i,j  ,k  ,n)+crse(i  ,j+1,k  ,n) &
                     &                         +crse(i,j  ,k+1,n)+crse(i  ,j+1,k+1,n))
                fine(ii+1,jj+1,kk+1,n) = 0.125d0*(crse(i,j  ,k  ,n)+crse(i+1,j  ,k  ,n) &
                     &                         +crse(i,j+1,k  ,n)+crse(i+1,j+1,k  ,n) &
                     &                         +crse(i,j  ,k+1,n)+crse(i+1,j  ,k+1,n) &
                     &                         +crse(i,j+1,k+1,n)+crse(i+1,j+1,k+1,n))
             end do
             i = chi(1)
             ii = i*2
             fine(ii  ,jj  ,kk  ,n) = crse(i,j,k,n)
             fine(ii  ,jj+1,kk  ,n) = 0.5d0  *(crse(i,j  ,k  ,n)+crse(i  ,j+1,k  ,n))
             fine(ii  ,jj  ,kk+1,n) = 0.5d0  *(crse(i,j  ,k  ,n)+crse(i  ,j  ,k+1,n))
             fine(ii  ,jj+1,kk+1,n) = 0.25d0 *(crse(i,j  ,k  ,n)+crse(i  ,j+1,k  ,n) &
                  &                         +crse(i,j  ,k+1,n)+crse(i  ,j+1,k+1,n))
          end do

          j = chi(2)
          jj = j*2
          do i = clo(1), chi(1)-1
             ii = i*2
             fine(ii  ,jj  ,kk  ,n) = crse(i,j,k,n)
             fine(ii+1,jj  ,kk  ,n) = 0.5d0  *(crse(i,j  ,k  ,n)+crse(i+1,j  ,k  ,n))
             fine(ii  ,jj  ,kk+1,n) = 0.5d0  *(crse(i,j  ,k  ,n)+crse(i  ,j  ,k+1,n))
             fine(ii+1,jj  ,kk+1,n) = 0.25d0 *(crse(i,j  ,k  ,n)+crse(i+1,j  ,k  ,n) &
                  &                         +crse(i,j  ,k+1,n)+crse(i+1,j  ,k+1,n))
          end do
          i = chi(1)
          ii = i*2
          fine(ii  ,jj  ,kk  ,n) = crse(i,j,k,n)
          fine(ii  ,jj  ,kk+1,n) = 0.5d0  *(crse(i,j  ,k  ,n)+crse(i  ,j  ,k+1,n))
       end do
    end do

    k = chi(3)
    kk = k*2

    do n = 1, nc
       do j = clo(2), chi(2)-1
          jj = j*2
          do i = clo(1), chi(1)-1
             ii = i*2
             fine(ii  ,jj  ,kk  ,n) = crse(i,j,k,n)
             fine(ii+1,jj  ,kk  ,n) = 0.5d0  *(crse(i,j  ,k  ,n)+crse(i+1,j  ,k  ,n))
             fine(ii  ,jj+1,kk  ,n) = 0.5d0  *(crse(i,j  ,k  ,n)+crse(i  ,j+1,k  ,n))
             fine(ii+1,jj+1,kk  ,n) = 0.25d0 *(crse(i,j  ,k  ,n)+crse(i+1,j  ,k  ,n) &
                  &                         +crse(i,j+1,k ,n )+crse(i+1,j+1,k ,n ))
          end do
          i = chi(1)
          ii = i*2
          fine(ii  ,jj  ,kk  ,n) = crse(i,j,k,n)
          fine(ii  ,jj+1,kk  ,n) = 0.5d0  *(crse(i,j  ,k ,n )+crse(i  ,j+1,k ,n ))
       end do
    end do
    
    j = chi(2)
    jj = j*2
    
    do n = 1, nc
       do i = clo(1), chi(1)-1
          ii = i*2
          fine(ii  ,jj  ,kk  ,n) = crse(i,j,k,n)
          fine(ii+1,jj  ,kk  ,n) = 0.5d0  *(crse(i,j  ,k  ,n)+crse(i+1,j  ,k  ,n))
       end do
       i = chi(1)
       ii = i*2
       fine(ii  ,jj  ,kk  ,n) = crse(i,j,k,n)
    end do

  end subroutine amrex_mlmg_lin_nd_interp

end module amrex_mlmg_interp_module
