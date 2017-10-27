
module amrex_mlpoisson_3d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlpoisson_adotx, amrex_mlpoisson_flux, amrex_mlpoisson_gsrb

contains

  subroutine amrex_mlpoisson_adotx (lo, hi, y, ylo, yhi, x, xlo, xhi, dxinv) &
       bind(c,name='amrex_mlpoisson_adotx')
    integer, dimension(3), intent(in) :: lo, hi, ylo, yhi, xlo, xhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    real(amrex_real), intent(in   ) :: x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    
    integer :: i,j,k
    real(amrex_real) :: dhx, dhy, dhz

    dhx = dxinv(1)*dxinv(1)
    dhy = dxinv(2)*dxinv(2)
    dhz = dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             y(i,j,k) = dhx * (x(i-1,j,k) - 2.d0*x(i,j,k) + x(i+1,j,k)) &
                  +     dhy * (x(i,j-1,k) - 2.d0*x(i,j,k) + x(i,j+1,k)) &
                  +     dhz * (x(i,j,k-1) - 2.d0*x(i,j,k) + x(i,j,k+1))
          end do
       end do
    end do
  end subroutine amrex_mlpoisson_adotx

  
  subroutine amrex_mlpoisson_flux (lo, hi, fx, fxlo, fxhi, fy, fylo, fyhi, &
       fz, fzlo, fzhi, sol, slo, shi, dxinv, face_only) &
       bind(c, name='amrex_mlpoisson_flux')
    integer, dimension(3), intent(in) :: lo, hi, fxlo, fxhi, fylo, fyhi, fzlo, fzhi, slo, shi
    real(amrex_real) :: dxinv(3)
    integer, value, intent(in) :: face_only
    real(amrex_real), intent(inout) :: fx (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
    real(amrex_real), intent(inout) :: fy (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
    real(amrex_real), intent(inout) :: fz (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
    real(amrex_real), intent(in   ) :: sol( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))    
    
    integer :: i,j,k
    real(amrex_real) :: dhx, dhy, dhz

    dhx = dxinv(1)
    dhy = dxinv(2)
    dhz = dxinv(3)

    if (face_only .eq. 1) then
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)+1, hi(1)+1-lo(1)
                fx(i,j,k) = dhx * (sol(i,j,k) - sol(i-1,j,k))
             end do
          end do
       end do
       
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)+1, hi(2)+1-lo(2)
             do i = lo(1), hi(1)
                fy(i,j,k) = dhy * (sol(i,j,k) - sol(i,j-1,k))
             end do
          end do
       end do
       
       do       k = lo(3), hi(3)+1, hi(3)+1-lo(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                fz(i,j,k) = dhz * (sol(i,j,k) - sol(i,j,k-1))
             end do
          end do
       end do

    else

       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)+1
                fx(i,j,k) = dhx * (sol(i,j,k) - sol(i-1,j,k))
             end do
          end do
       end do
       
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)+1
             do i = lo(1), hi(1)
                fy(i,j,k) = dhy * (sol(i,j,k) - sol(i,j-1,k))
             end do
          end do
       end do
       
       do       k = lo(3), hi(3)+1
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                fz(i,j,k) = dhz * (sol(i,j,k) - sol(i,j,k-1))
             end do
          end do
       end do
    end if

  end subroutine amrex_mlpoisson_flux


  subroutine amrex_mlpoisson_gsrb (lo, hi, phi, hlo, hhi, rhs, rlo, rhi, &
       f0,f0lo,f0hi,f1,f1lo,f1hi,f2,f2lo,f2hi,f3,f3lo,f3hi,f4,f4lo,f4hi,f5,f5lo,f5hi, &
       m0,m0lo,m0hi,m1,m1lo,m1hi,m2,m2lo,m2hi,m3,m3lo,m3hi,m4,m4lo,m4hi,m5,m5lo,m5hi, &
       blo, bhi, dxinv, redblack) bind(c,name='amrex_mlpoisson_gsrb')
    integer, dimension(3), intent(in) :: lo, hi, hlo, hhi, rlo, rhi, &
         f0lo,f0hi,f1lo,f1hi,f2lo,f2hi,f3lo,f3hi,f4lo,f4hi,f5lo,f5hi, &
         m0lo,m0hi,m1lo,m1hi,m2lo,m2hi,m3lo,m3hi,m4lo,m4hi,m5lo,m5hi, blo, bhi
    integer, intent(in), value :: redblack
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: phi(hlo(1):hhi(1),hlo(2):hhi(2),hlo(3):hhi(3))
    real(amrex_real), intent(in   ) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in   ) :: f0(f0lo(1):f0hi(1),f0lo(2):f0hi(2),f0lo(3):f0hi(3))
    real(amrex_real), intent(in   ) :: f1(f1lo(1):f1hi(1),f1lo(2):f1hi(2),f1lo(3):f1hi(3))
    real(amrex_real), intent(in   ) :: f2(f2lo(1):f2hi(1),f2lo(2):f2hi(2),f2lo(3):f2hi(3))
    real(amrex_real), intent(in   ) :: f3(f3lo(1):f3hi(1),f3lo(2):f3hi(2),f3lo(3):f3hi(3))
    real(amrex_real), intent(in   ) :: f4(f4lo(1):f4hi(1),f4lo(2):f4hi(2),f4lo(3):f4hi(3))
    real(amrex_real), intent(in   ) :: f5(f5lo(1):f5hi(1),f5lo(2):f5hi(2),f5lo(3):f5hi(3))
    integer         , intent(in   ) :: m0(m0lo(1):m0hi(1),m0lo(2):m0hi(2),m0lo(3):m0hi(3))
    integer         , intent(in   ) :: m1(m1lo(1):m1hi(1),m1lo(2):m1hi(2),m1lo(3):m1hi(3))
    integer         , intent(in   ) :: m2(m2lo(1):m2hi(1),m2lo(2):m2hi(2),m2lo(3):m2hi(3))
    integer         , intent(in   ) :: m3(m3lo(1):m3hi(1),m3lo(2):m3hi(2),m3lo(3):m3hi(3))
    integer         , intent(in   ) :: m4(m4lo(1):m4hi(1),m4lo(2):m4hi(2),m4lo(3):m4hi(3))
    integer         , intent(in   ) :: m5(m5lo(1):m5hi(1),m5lo(2):m5hi(2),m5lo(3):m5hi(3))

    integer :: i,j, k, ioff
    real(amrex_real) :: dhx, dhy, dhz, cf0, cf1, cf2, cf3, cf4, cf5
    real(amrex_real) :: gamma, g_m_d, res
    real(amrex_real), parameter :: omega = 1.15d0

    dhx = dxinv(1)*dxinv(1)
    dhy = dxinv(2)*dxinv(2)
    dhz = dxinv(3)*dxinv(3)

    gamma = -2.d0*(dhx+dhy+dhz)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          ioff = mod(lo(1) + j + k + redblack,2)
          do i = lo(1) + ioff, hi(1), 2

             cf0 = merge(f0(blo(1),j,k), 0.d0,  &
                  &      (i .eq. blo(1)) .and. (m0(blo(1)-1,j,k).gt.0))
             cf1 = merge(f1(i,blo(2),k), 0.d0,  &
                  &      (j .eq. blo(2)) .and. (m1(i,blo(2)-1,k).gt.0))
             cf2 = merge(f2(i,j,blo(3)), 0.d0,  &
                  &      (k .eq. blo(3)) .and. (m2(i,j,blo(3)-1).gt.0))
             cf3 = merge(f3(bhi(1),j,k), 0.d0,  &
                  &      (i .eq. bhi(1)) .and. (m3(bhi(1)+1,j,k).gt.0))
             cf4 = merge(f4(i,bhi(2),k), 0.d0,  &
                  &      (j .eq. bhi(2)) .and. (m4(i,bhi(2)+1,k).gt.0))
             cf5 = merge(f5(i,j,bhi(3)), 0.d0,  &
                  &      (k .eq. bhi(3)) .and. (m5(i,j,bhi(3)+1).gt.0))

             g_m_d = gamma + dhx*(cf0+cf3) + dhy*(cf1+cf4) + dhz*(cf2+cf5)

             res = rhs(i,j,k) - gamma*phi(i,j,k) &
                  - dhx*(phi(i-1,j,k) + phi(i+1,j,k))  &
                  - dhy*(phi(i,j-1,k) + phi(i,j+1,k))  &
                  - dhz*(phi(i,j,k-1) + phi(i,j,k+1))

             phi(i,j,k) = phi(i,j,k) + omega/g_m_d * res

          end do
       end do
    end do

  end subroutine amrex_mlpoisson_gsrb

end module amrex_mlpoisson_3d_module
