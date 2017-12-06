module amrex_mlnodelap_3d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlndlap_avgdown_coeff, amrex_mlndlap_divu, &
       amrex_mlndlap_adotx, amrex_mlndlap_jacobi

contains

  subroutine amrex_mlndlap_avgdown_coeff (lo, hi, crse, clo, chi, fine, flo, fhi, idim) &
       bind(c,name='amrex_mlndlap_avgdown_coeff')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, flo, fhi
    integer, intent(in) :: idim
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
  end subroutine amrex_mlndlap_avgdown_coeff


  subroutine amrex_mlndlap_divu (lo, hi, rhs, rlo, rhi, vel, vlo, vhi, dxinv) &
       bind(c,name='amrex_mlndlap_divu')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, vlo, vhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)
  end subroutine amrex_mlndlap_divu


  subroutine amrex_mlndlap_adotx (lo, hi, y, ylo, yhi, x, xlo, xhi, &
       sx, sxlo, sxhi, sy, sylo, syhi, dxinv) bind(c,name='amrex_mlndlap_adotx')
    integer, dimension(2), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, sxlo, sxhi, sylo, syhi
    real(amrex_real), intent(in) :: dxinv(2)
    real(amrex_real), intent(inout) ::  y( ylo(1): yhi(1), ylo(2): yhi(2))
    real(amrex_real), intent(in   ) ::  x( xlo(1): xhi(1), xlo(2): xhi(2))
    real(amrex_real), intent(in   ) :: sx(sxlo(1):sxhi(1),sxlo(2):sxhi(2))
    real(amrex_real), intent(in   ) :: sy(sylo(1):syhi(1),sylo(2):syhi(2))
    
    integer :: i,j
!    real(amrex_real) :: dhx, dhy

  end subroutine amrex_mlndlap_adotx  


  subroutine amrex_mlndlap_jacobi (lo, hi, sol, slo, shi, Ax, alo, ahi, rhs, rlo, rhi, &
       sx, sxlo, sxhi, sy, sylo, syhi, sz, szlo, szhi, dxinv) bind(c,name='amrex_mlndlap_jacobi')
    integer, dimension(3),intent(in) :: lo,hi,slo,shi,alo,ahi,rlo,rhi,sxlo,sxhi,sylo,syhi,szlo,szhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    real(amrex_real), intent(in   ) :: Ax ( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: sx (sxlo(1):sxhi(1),sxlo(2):sxhi(2),sxlo(3):sxhi(3))
    real(amrex_real), intent(in   ) :: sy (sylo(1):syhi(1),sylo(2):syhi(2),sylo(3):syhi(3))
    real(amrex_real), intent(in   ) :: sz (szlo(1):szhi(1),szlo(2):szhi(2),szlo(3):szhi(3))
    
  end subroutine amrex_mlndlap_jacobi

end module amrex_mlnodelap_3d_module
