
module mllinop_2d_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mllinop_grad

contains
  
  subroutine amrex_mllinop_grad (xlo, xhi, ylo, yhi, dat, dlo, dhi, gx, gxlo, gxhi, &
       gy, gylo, gyhi, dxinv) bind(c,name='amrex_mllinop_grad')
    integer, dimension(2), intent(in) :: xlo, xhi, ylo, yhi, dlo, dhi, gxlo, gxhi, gylo, gyhi
    real(amrex_real), dimension(2), intent(in) :: dxinv
    real(amrex_real), intent(in   ) :: dat( dlo(1): dhi(1), dlo(2): dhi(2))
    real(amrex_real), intent(inout) :: gx (gxlo(1):gxhi(1),gxlo(2):gxhi(2))
    real(amrex_real), intent(inout) :: gy (gylo(1):gyhi(1),gylo(2):gyhi(2))
    
    integer :: i,j
    real(amrex_real) :: dhx, dhy
    
    dhx = dxinv(1)
    dhy = dxinv(2)

    do    j = xlo(2), xhi(2)
       do i = xlo(1), xhi(1)
          gx(i,j) = dhx * (dat(i,j) - dat(i-1,j))
       end do
    end do
    
    do    j = ylo(2), yhi(2)
       do i = ylo(1), yhi(1)
          gy(i,j) = dhy * (dat(i,j) - dat(i,j-1))
       end do
    end do

  end subroutine amrex_mllinop_grad

end module mllinop_2d_module
