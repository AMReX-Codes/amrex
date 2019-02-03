
module amrex_mlabeclap_1d_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mlabeclap_flux

contains

  subroutine amrex_mlabeclap_flux (lo, hi, fx, fxlo, fxhi, sol, slo, shi, bx, bxlo, bxhi, &
       dxinv, beta, face_only) bind(c, name='amrex_mlabeclap_flux')
    integer, dimension(1), intent(in) :: lo, hi, fxlo, fxhi, slo, shi, bxlo, bxhi
    real(amrex_real) :: dxinv(1)
    real(amrex_real), value, intent(in) :: beta
    integer, value, intent(in) :: face_only
    real(amrex_real), intent(inout) :: fx (fxlo(1):fxhi(1))
    real(amrex_real), intent(in   ) :: sol( slo(1): shi(1))
    real(amrex_real), intent(in   ) :: bx (bxlo(1):bxhi(1))
    
    integer :: i
    real(amrex_real) :: dhx

    dhx = beta*dxinv(1)

    if (face_only .eq. 1) then
       do i = lo(1), hi(1)+1, hi(1)+1-lo(1)
          fx(i) = -dhx * bx(i)*(sol(i) - sol(i-1))
       end do       
    else
       do i = lo(1), hi(1)+1
          fx(i) = -dhx * bx(i)*(sol(i) - sol(i-1))
       end do
    end if
  end subroutine amrex_mlabeclap_flux

end module amrex_mlabeclap_1d_module
