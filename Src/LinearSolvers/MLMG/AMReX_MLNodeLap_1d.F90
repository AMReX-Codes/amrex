module amrex_mlnodelap_1d_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use amrex_constants_module

  implicit none

  private
  public :: &
       ! masks
       amrex_mlndlap_any_fine_sync_cells, &
       ! coeffs
       ! bc
       ! operator
       ! restriction
       ! interpolation
       ! rhs & u
       ! residual
       amrex_mlndlap_crse_resid, &
       amrex_mlndlap_res_fine_contrib, amrex_mlndlap_res_cf_contrib, &
       ! sync residual
       amrex_mlndlap_zero_fine

  ! RAP

contains

  function amrex_mlndlap_any_fine_sync_cells (lo, hi, msk, mlo, mhi, fine_flag) result(r) &
       bind(c,name='amrex_mlndlap_any_fine_sync_cells')
    integer :: r
    integer, dimension(1), intent(in) :: lo, hi, mlo, mhi
    integer, intent(in   ) :: msk  ( mlo(1): mhi(1))
    integer, intent(in) :: fine_flag
  end function amrex_mlndlap_any_fine_sync_cells


  subroutine amrex_mlndlap_crse_resid (lo, hi, resid, rslo, rshi, rhs, rhlo, rhhi, msk, mlo, mhi, &
       ndlo, ndhi, bclo, bchi) bind(c, name='amrex_mlndlap_crse_resid')
    integer, dimension(1), intent(in) :: lo, hi, rslo, rshi, rhlo, rhhi, mlo, mhi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(inout) :: resid(rslo(1):rshi(1))
    real(amrex_real), intent(in   ) :: rhs  (rhlo(1):rhhi(1))
    integer         , intent(in   ) :: msk  ( mlo(1): mhi(1))
  end subroutine amrex_mlndlap_crse_resid


  subroutine amrex_mlndlap_res_fine_contrib (clo, chi, cglo, cghi, f, flo, fhi, &
       x, xlo, xhi, sig, slo, shi, Ax, alo, ahi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_res_fine_contrib')
    integer, dimension(1), intent(in) :: clo, chi, cglo, cghi, flo, fhi, xlo, xhi, &
         slo, shi, alo, ahi, mlo, mhi
    real(amrex_real), intent(inout) :: f  (flo(1):fhi(1))
    real(amrex_real), intent(in   ) :: x  (xlo(1):xhi(1))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1))
    real(amrex_real), intent(inout) :: Ax (alo(1):ahi(1))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1))
    real(amrex_real), intent(in) :: dxinv(1)
  end subroutine amrex_mlndlap_res_fine_contrib


  subroutine amrex_mlndlap_res_cf_contrib (lo, hi, res, rlo, rhi, phi, phlo, phhi, &
       rhs, rhlo, rhhi, sig, slo, shi, dmsk, mlo, mhi, ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi, &
       fc, clo, chi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_res_cf_contrib')
    integer, dimension(1), intent(in) :: lo, hi, rlo, rhi, phlo, phhi, rhlo, rhhi, slo, shi, &
         mlo, mhi, nmlo, nmhi, cmlo, cmhi, clo, chi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: res( rlo(1): rhi(1))
    real(amrex_real), intent(in   ) :: phi(phlo(1):phhi(1))
    real(amrex_real), intent(in   ) :: rhs(rhlo(1):rhhi(1))
    real(amrex_real), intent(in   ) :: sig( slo(1): shi(1))
    real(amrex_real), intent(inout) :: fc ( clo(1): chi(1))
    integer, intent(in) :: dmsk(mlo(1):mhi(1))
    integer, intent(in) :: ndmsk(nmlo(1):nmhi(1))
    integer, intent(in) :: ccmsk(cmlo(1):cmhi(1))
  end subroutine amrex_mlndlap_res_cf_contrib


  subroutine amrex_mlndlap_zero_fine (lo, hi, phi, dlo, dhi, msk, mlo, mhi, fine_flag) &
       bind(c, name='amrex_mlndlap_zero_fine')
    integer, dimension(1), intent(in) :: lo, hi, dlo, dhi, mlo, mhi
    real(amrex_real), intent(inout) :: phi(dlo(1):dhi(1))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1))
    integer, intent(in) :: fine_flag
  end subroutine amrex_mlndlap_zero_fine

end module amrex_mlnodelap_1d_module
