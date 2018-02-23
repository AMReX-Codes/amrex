
module amrex_mlpoisson_2d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlpoisson_adotx, amrex_mlpoisson_normalize, amrex_mlpoisson_flux, amrex_mlpoisson_gsrb

contains

  subroutine amrex_mlpoisson_adotx (lo, hi, y, ylo, yhi, x, xlo, xhi, rc, re, rlo, rhi, dxinv) &
       bind(c,name='amrex_mlpoisson_adotx')
    integer, dimension(2), intent(in) :: lo, hi, ylo, yhi, xlo, xhi
    integer, intent(in) :: rlo, rhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: y(ylo(1):yhi(1),ylo(2):yhi(2))
    real(amrex_real), intent(in   ) :: x(xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in) :: rc(rlo:rhi)
    real(amrex_real), intent(in) :: re(rlo:rhi+1)
    
    integer :: i,j
    real(amrex_real) :: dhx, dhy

    dhx = dxinv(1)*dxinv(1)
    dhy = dxinv(2)*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          y(i,j) = dhx * (re(i)*x(i-1,j) - (re(i)+re(i+1))*x(i,j) + re(i+1)*x(i+1,j)) &
               +   dhy * rc(i) * (x(i,j-1) - 2.d0*x(i,j) + x(i,j+1))
       end do
    end do
  end subroutine amrex_mlpoisson_adotx


  subroutine amrex_mlpoisson_normalize (lo, hi, x, xlo, xhi, rc, re, rlo, rhi, dxinv) &
       bind(c,name='amrex_mlpoisson_normalize')
    integer, dimension(2), intent(in) :: lo, hi, xlo, xhi
    integer, intent(in) :: rlo, rhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: x(xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in) :: rc(rlo:rhi)
    real(amrex_real), intent(in) :: re(rlo:rhi+1)
    
    integer :: i,j
    real(amrex_real) :: dhx, dhy

    dhx = dxinv(1)*dxinv(1)
    dhy = dxinv(2)*dxinv(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          x(i,j) = x(i,j) / (-dhx*(re(i)+re(i+1)) - dhy*rc(i)*2.d0)
       end do
    end do
  end subroutine amrex_mlpoisson_normalize

  
  subroutine amrex_mlpoisson_flux (lo, hi, fx, fxlo, fxhi, fy, fylo, fyhi, &
       sol, slo, shi, rc, re, rlo, rhi, dxinv, face_only) bind(c, name='amrex_mlpoisson_flux')
    integer, dimension(2), intent(in) :: lo, hi, fxlo, fxhi, fylo, fyhi, slo, shi
    integer, intent(in) :: rlo, rhi
    real(amrex_real) :: dxinv(2)
    integer, value, intent(in) :: face_only
    real(amrex_real), intent(inout) :: fx (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(amrex_real), intent(inout) :: fy (fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(amrex_real), intent(in   ) :: sol( slo(1): shi(1), slo(2): shi(2))
    real(amrex_real), intent(in) :: rc(rlo:rhi)
    real(amrex_real), intent(in) :: re(rlo:rhi+1)
    
    integer :: i,j
    real(amrex_real) :: dhx, dhy

    dhx = dxinv(1)
    dhy = dxinv(2)

    if (face_only .eq. 1) then
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)+1, hi(1)+1-lo(1)
             fx(i,j) = dhx * re(i) * (sol(i,j) - sol(i-1,j))
          end do
       end do
       
       do    j = lo(2), hi(2)+1, hi(2)+1-lo(2)
          do i = lo(1), hi(1)
             fy(i,j) = dhy * rc(i) * (sol(i,j) - sol(i,j-1))
          end do
       end do

    else

       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             fx(i,j) = dhx * re(i) * (sol(i,j) - sol(i-1,j))
          end do
       end do
       
       do    j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             fy(i,j) = dhy * rc(i) * (sol(i,j) - sol(i,j-1))
          end do
       end do
    end if

  end subroutine amrex_mlpoisson_flux


  subroutine amrex_mlpoisson_gsrb (lo, hi, phi, hlo, hhi, rhs, rlo, rhi, &
       f0,f0lo,f0hi,f1,f1lo,f1hi,f2,f2lo,f2hi,f3,f3lo,f3hi, &
       m0,m0lo,m0hi,m1,m1lo,m1hi,m2,m2lo,m2hi,m3,m3lo,m3hi, &
       rc, re, blo, bhi, dxinv, redblack) bind(c,name='amrex_mlpoisson_gsrb')
    integer, dimension(2), intent(in) :: lo, hi, hlo, hhi, rlo, rhi, &
         f0lo,f0hi,f1lo,f1hi,f2lo,f2hi,f3lo,f3hi, &
         m0lo,m0hi,m1lo,m1hi,m2lo,m2hi,m3lo,m3hi, blo, bhi
    integer, intent(in), value :: redblack
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: phi(hlo(1):hhi(1),hlo(2):hhi(2))
    real(amrex_real), intent(in   ) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) :: f0(f0lo(1):f0hi(1),f0lo(2):f0hi(2))
    real(amrex_real), intent(in   ) :: f1(f1lo(1):f1hi(1),f1lo(2):f1hi(2))
    real(amrex_real), intent(in   ) :: f2(f2lo(1):f2hi(1),f2lo(2):f2hi(2))
    real(amrex_real), intent(in   ) :: f3(f3lo(1):f3hi(1),f3lo(2):f3hi(2))
    integer         , intent(in   ) :: m0(m0lo(1):m0hi(1),m0lo(2):m0hi(2))
    integer         , intent(in   ) :: m1(m1lo(1):m1hi(1),m1lo(2):m1hi(2))
    integer         , intent(in   ) :: m2(m2lo(1):m2hi(1),m2lo(2):m2hi(2))
    integer         , intent(in   ) :: m3(m3lo(1):m3hi(1),m3lo(2):m3hi(2))
    real(amrex_real), intent(in) :: rc(blo(1):bhi(1))
    real(amrex_real), intent(in) :: re(blo(1):bhi(1)+1)

    integer :: i,j, ioff
    real(amrex_real) :: dhx, dhy, cf0, cf1, cf2, cf3
    real(amrex_real) :: gamma, g_m_d, res

    dhx = dxinv(1)*dxinv(1)
    dhy = dxinv(2)*dxinv(2)

    do j = lo(2), hi(2)
       ioff = mod(lo(1) + j + redblack,2)
       do i = lo(1) + ioff, hi(1), 2

          cf0 = merge(f0(blo(1),j), 0.d0,  &
               &      (i .eq. blo(1)) .and. (m0(blo(1)-1,j).gt.0))
          cf1 = merge(f1(i,blo(2)), 0.d0,  &
               &      (j .eq. blo(2)) .and. (m1(i,blo(2)-1).gt.0))
          cf2 = merge(f2(bhi(1),j), 0.d0,  &
               &      (i .eq. bhi(1)) .and. (m2(bhi(1)+1,j).gt.0))
          cf3 = merge(f3(i,bhi(2)), 0.d0,  &
               &      (j .eq. bhi(2)) .and. (m3(i,bhi(2)+1).gt.0))

          gamma = -dhx*(re(i)+re(i+1)) - 2.d0*dhy*rc(i)

          g_m_d = gamma + dhx*(re(i)*cf0+re(i+1)*cf2) + dhy*rc(i)*(cf1+cf3)

          res = rhs(i,j) - gamma*phi(i,j) &
               - dhx*(re(i)*phi(i-1,j) + re(i+1)*phi(i+1,j))  &
               - dhy*rc(i)*(phi(i,j-1) + phi(i,j+1))

          phi(i,j) = phi(i,j) + res /g_m_d

       end do
    end do

  end subroutine amrex_mlpoisson_gsrb

end module amrex_mlpoisson_2d_module
