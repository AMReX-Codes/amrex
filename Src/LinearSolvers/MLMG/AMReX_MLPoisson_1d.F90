
module amrex_mlpoisson_1d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlpoisson_adotx, amrex_mlpoisson_normalize, amrex_mlpoisson_flux, amrex_mlpoisson_gsrb

contains

  subroutine amrex_mlpoisson_adotx (lo, hi, y, ylo, yhi, x, xlo, xhi, rc, re, rlo, rhi, dxinv) &
       bind(c,name='amrex_mlpoisson_adotx')
    integer, dimension(1), intent(in) :: lo, hi, ylo, yhi, xlo, xhi
    integer, intent(in) :: rlo, rhi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: y(ylo(1):yhi(1))
    real(amrex_real), intent(in   ) :: x(xlo(1):xhi(1))
    real(amrex_real), intent(in) :: rc(rlo:rhi)
    real(amrex_real), intent(in) :: re(rlo:rhi+1)
    
    integer :: i
    real(amrex_real) :: dhx

    dhx = dxinv(1)*dxinv(1)

    do i = lo(1), hi(1)
       y(i) = dhx * (re(i)*x(i-1) - (re(i)+re(i+1))*x(i) + re(i+1)*x(i+1))
    end do
  end subroutine amrex_mlpoisson_adotx


  subroutine amrex_mlpoisson_normalize (lo, hi, x, xlo, xhi, rc, re, rlo, rhi, dxinv) &
       bind(c,name='amrex_mlpoisson_normalize')
    integer, dimension(1), intent(in) :: lo, hi, xlo, xhi
    integer, intent(in) :: rlo, rhi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: x(xlo(1):xhi(1))
    real(amrex_real), intent(in) :: rc(rlo:rhi)
    real(amrex_real), intent(in) :: re(rlo:rhi+1)
    
    integer :: i
    real(amrex_real) :: dhx

    dhx = dxinv(1)*dxinv(1)

    do i = lo(1), hi(1)
       x(i) = x(i) / (-dhx*(re(i)+re(i+1)))
    end do
  end subroutine amrex_mlpoisson_normalize

  
  subroutine amrex_mlpoisson_flux (lo, hi, fx, fxlo, fxhi, sol, slo, shi, rc, re, rlo, rhi, &
       dxinv, face_only) bind(c, name='amrex_mlpoisson_flux')
    integer, dimension(1), intent(in) :: lo, hi, fxlo, fxhi, slo, shi
    integer, intent(in) :: rlo, rhi
    real(amrex_real) :: dxinv(1)
    integer, value, intent(in) :: face_only
    real(amrex_real), intent(inout) :: fx (fxlo(1):fxhi(1))
    real(amrex_real), intent(in   ) :: sol( slo(1): shi(1))
    real(amrex_real), intent(in) :: rc(rlo:rhi)
    real(amrex_real), intent(in) :: re(rlo:rhi+1)
    
    integer :: i
    real(amrex_real) :: dhx

    dhx = dxinv(1)

    if (face_only .eq. 1) then
       do i = lo(1), hi(1)+1, hi(1)+1-lo(1)
          fx(i) = dhx * re(i) * (sol(i) - sol(i-1))
       end do
    else
       do i = lo(1), hi(1)+1
          fx(i) = dhx * re(i) * (sol(i) - sol(i-1))
       end do
    end if

  end subroutine amrex_mlpoisson_flux


  subroutine amrex_mlpoisson_gsrb (lo, hi, phi, hlo, hhi, rhs, rlo, rhi, &
       f0,f0lo,f0hi,f1,f1lo,f1hi, m0,m0lo,m0hi,m1,m1lo,m1hi, &
       rc, re, blo, bhi, dxinv, redblack) bind(c,name='amrex_mlpoisson_gsrb')
    integer, dimension(1), intent(in) :: lo, hi, hlo, hhi, rlo, rhi, &
         f0lo,f0hi,f1lo,f1hi, m0lo,m0hi,m1lo,m1hi, blo, bhi
    integer, intent(in), value :: redblack
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: phi(hlo(1):hhi(1))
    real(amrex_real), intent(in   ) :: rhs(rlo(1):rhi(1))
    real(amrex_real), intent(in   ) :: f0(f0lo(1):f0hi(1))
    real(amrex_real), intent(in   ) :: f1(f1lo(1):f1hi(1))
    integer         , intent(in   ) :: m0(m0lo(1):m0hi(1))
    integer         , intent(in   ) :: m1(m1lo(1):m1hi(1))
    real(amrex_real), intent(in) :: rc(blo(1):bhi(1))
    real(amrex_real), intent(in) :: re(blo(1):bhi(1)+1)

    integer :: i, ioff
    real(amrex_real) :: dhx, cf0, cf1
    real(amrex_real) :: gamma, g_m_d, res

    dhx = dxinv(1)*dxinv(1)

    ioff = mod(lo(1) + redblack,2)
    do i = lo(1) + ioff, hi(1), 2

       cf0 = merge(f0(blo(1)), 0.d0,  &
            &      (i .eq. blo(1)) .and. (m0(blo(1)-1).gt.0))
       cf1 = merge(f1(bhi(1)), 0.d0,  &
            &      (i .eq. bhi(1)) .and. (m1(bhi(1)+1).gt.0))
       
       gamma = -dhx*(re(i)+re(i+1))
       
       g_m_d = gamma + dhx*(re(i)*cf0+re(i+1)*cf1)

       res = rhs(i) - gamma*phi(i) &
            - dhx*(re(i)*phi(i-1) + re(i+1)*phi(i+1))
       
       phi(i) = phi(i) + res /g_m_d

    end do

  end subroutine amrex_mlpoisson_gsrb

end module amrex_mlpoisson_1d_module
