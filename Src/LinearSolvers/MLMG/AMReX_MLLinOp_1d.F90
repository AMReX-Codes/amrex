
module mllinop_1d_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mllinop_grad

contains

  subroutine amrex_mllinop_grad (xlo, xhi, dat, dlo, dhi, gx, gxlo, gxhi, dxinv) &
       bind(c,name='amrex_mllinop_grad')
    integer, dimension(1), intent(in) :: xlo, xhi, dlo, dhi, gxlo, gxhi
    real(amrex_real), dimension(1), intent(in) :: dxinv
    real(amrex_real), intent(in   ) :: dat( dlo(1): dhi(1))
    real(amrex_real), intent(inout) :: gx (gxlo(1):gxhi(1))
    
    integer :: i
    real(amrex_real) :: dhx

    dhx = dxinv(1)

    do i = xlo(1), xhi(1)
       gx(i) = dhx * (dat(i) - dat(i-1))
    end do
  end subroutine amrex_mllinop_grad

end module mllinop_1d_module
