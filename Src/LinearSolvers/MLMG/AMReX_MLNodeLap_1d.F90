module amrex_mlnodelap_1d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlndlap_sigma_cctoedge, amrex_mlndlap_divu

contains

  subroutine amrex_mlndlap_sigma_cctoedge (xlo, xhi, sigx, sxlo, sxhi, sigcc, clo, chi) &
       bind(c, name='amrex_mlndlap_sigma_cctoedge')
    integer, dimension(1), intent(in) :: xlo, xhi, sxlo, sxhi, clo, chi
    real(amrex_real), intent(inout) :: sigx (sxlo(1):sxhi(1))
    real(amrex_real), intent(in   ) :: sigcc( clo(1): chi(1))
    integer :: i
    do i = xlo(1), xhi(1)
       sigx(i) = sigcc(i)
    end do
  end subroutine amrex_mlndlap_sigma_cctoedge


  subroutine amrex_mlndlap_divu (lo, hi, rhs, rlo, rhi, vel, vlo, vhi, dxinv) &
       bind(c,name='amrex_mlndlap_divu')
    integer, dimension(1), intent(in) :: lo, hi, rlo, rhi, vlo, vhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1))
  end subroutine amrex_mlndlap_divu

end module amrex_mlnodelap_1d_module
