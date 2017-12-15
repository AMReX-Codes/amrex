
module mllinop_3d_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mllinop_grad

contains

  subroutine amrex_mllinop_grad (xlo, xhi, ylo, yhi, zlo, zhi, dat, dlo, dhi, gx, gxlo, gxhi, &
       gy, gylo, gyhi, gz, gzlo, gzhi, dxinv) bind(c,name='amrex_mllinop_grad')
    integer, dimension(3), intent(in) :: xlo, xhi, ylo, yhi, zlo, zhi, dlo, dhi, &
         gxlo, gxhi, gylo, gyhi, gzlo, gzhi
    real(amrex_real), dimension(3), intent(in) :: dxinv
    real(amrex_real), intent(in   ) :: dat( dlo(1): dhi(1), dlo(2): dhi(2), dlo(3): dhi(3))
    real(amrex_real), intent(inout) :: gx (gxlo(1):gxhi(1),gxlo(2):gxhi(2),gxlo(3):gxhi(3))
    real(amrex_real), intent(inout) :: gy (gylo(1):gyhi(1),gylo(2):gyhi(2),gylo(3):gyhi(3))
    real(amrex_real), intent(inout) :: gz (gzlo(1):gzhi(1),gzlo(2):gzhi(2),gzlo(3):gzhi(3))

    integer :: i,j, k
    real(amrex_real) :: dhx, dhy, dhz
    
    dhx = dxinv(1)
    dhy = dxinv(2)
    dhz = dxinv(3)

    do       k = xlo(3), xhi(3)
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             gx(i,j,k) = dhx * (dat(i,j,k) - dat(i-1,j,k))
          end do
       end do
    end do
    
    do       k = ylo(3), yhi(3)
       do    j = ylo(2), yhi(2)
          do i = ylo(1), yhi(1)
             gy(i,j,k) = dhy * (dat(i,j,k) - dat(i,j-1,k))
          end do
       end do
    end do

    do       k = zlo(3), zhi(3)
       do    j = zlo(2), zhi(2)
          do i = zlo(1), zhi(1)
             gz(i,j,k) = dhz * (dat(i,j,k) - dat(i,j,k-1))
          end do
       end do
    end do

  end subroutine amrex_mllinop_grad

end module mllinop_3d_module
