
module amrex_mlpoisson_3d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlpoisson_adotx, amrex_mlpoisson_flux

contains

  subroutine amrex_mlpoisson_adotx (lo, hi, y, ylo, yhi, x, xlo, xhi, dxinv) &
       bind(c,name='amrex_mlpoisson_adotx')
    integer, dimension(3), intent(in) :: lo, hi, ylo, yhi, xlo, xhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) ::  y( ylo(1): yhi(1), ylo(2): yhi(2), ylo(3): yhi(3))
    real(amrex_real), intent(in   ) ::  x( xlo(1): xhi(1), xlo(2): xhi(2), xlo(3): xhi(3))
    
    integer :: i,j,k
    real(amrex_real) :: dhx, dhy, dhz

    dhx = dxinv(1)*dxinv(1)
    dhy = dxinv(2)*dxinv(2)
    dhz = dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             y(i,j,k) = dhx * (x(i-1,j,k) + 2.d0*x(i,j,k) + x(i+1,j,k)) &
                  +     dhy * (x(i,j-1,k) + 2.d0*x(i,j,k) + x(i,j+1,k)) &
                  +     dhz * (x(i,j,k-1) + 2.d0*x(i,j,k) + x(i,j,k+1))
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

end module amrex_mlpoisson_3d_module
