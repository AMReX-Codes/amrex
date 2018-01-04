module amrex_mlnodelap_1d_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann, amrex_lo_inflow, amrex_lo_periodic
  implicit none

  ! external dirichlet at physical boundary or internal dirichlet at crse/fine boundary
  integer, parameter :: dirichlet = 1

  integer, parameter :: crse_cell = 0
  integer, parameter :: fine_cell = 1
  integer, parameter :: crse_node = 0
  integer, parameter :: crse_fine_node = 1
  integer, parameter :: fine_node = 2

  private
  public :: &
       ! masks
       amrex_mlndlap_set_nodal_mask, amrex_mlndlap_set_dirichlet_mask, &
       amrex_mlndlap_fixup_res_mask, amrex_mlndlap_set_dot_mask, &
       amrex_mlndlap_any_fine_sync_cells, &
       ! coeffs
       amrex_mlndlap_avgdown_coeff, amrex_mlndlap_fillbc_cc, &
       ! bc
       amrex_mlndlap_applybc, &
       ! operator
       amrex_mlndlap_adotx_ha, amrex_mlndlap_adotx_aa, &
       amrex_mlndlap_jacobi_ha, amrex_mlndlap_jacobi_aa, &
       amrex_mlndlap_gauss_seidel_ha, amrex_mlndlap_gauss_seidel_aa, &
       ! restriction
       amrex_mlndlap_restriction, &
       ! interpolation
       amrex_mlndlap_interpolation_ha, amrex_mlndlap_interpolation_aa, &
       ! rhs & u
       amrex_mlndlap_divu, amrex_mlndlap_add_rhcc, amrex_mlndlap_mknewu, &
       amrex_mlndlap_divu_fine_contrib, amrex_mlndlap_divu_cf_contrib, &
       ! residual
       amrex_mlndlap_crse_resid, &
       amrex_mlndlap_res_fine_contrib, amrex_mlndlap_res_cf_contrib, &
       ! sync residual
       amrex_mlndlap_zero_fine

contains

  subroutine amrex_mlndlap_set_nodal_mask (lo, hi, nmsk, nlo, nhi, cmsk, clo, chi) &
       bind(c,name='amrex_mlndlap_set_nodal_mask')
    integer, dimension(1), intent(in) :: lo, hi, nlo, nhi, clo, chi
    integer, intent(inout) :: nmsk(nlo(1):nhi(1))
    integer, intent(in   ) :: cmsk(clo(1):chi(1))
  end subroutine amrex_mlndlap_set_nodal_mask


  subroutine amrex_mlndlap_set_dirichlet_mask (dmsk, dlo, dhi, omsk, olo, ohi, &
       domlo, domhi, bclo, bchi) bind(c,name='amrex_mlndlap_set_dirichlet_mask')
    integer, dimension(1) :: dlo, dhi, olo, ohi, domlo, domhi, bclo, bchi
    integer, intent(inout) :: dmsk(dlo(1):dhi(1))
    integer, intent(in   ) :: omsk(olo(1):ohi(1))
  end subroutine amrex_mlndlap_set_dirichlet_mask


  subroutine amrex_mlndlap_fixup_res_mask (lo, hi, rmsk, rlo, rhi, fmsk, flo, fhi) &
       bind(c,name='amrex_mlndlap_fixup_res_mask')
    integer, dimension(1), intent(in) :: lo, hi, rlo, rhi, flo, fhi
    integer, intent(inout) :: rmsk(rlo(1):rhi(1))
    integer, intent(in   ) :: fmsk(flo(1):fhi(1))
  end subroutine amrex_mlndlap_fixup_res_mask


  subroutine amrex_mlndlap_set_dot_mask (lo, hi, dmsk, dlo, dhi, omsk, olo, ohi, &
       domlo, domhi, bclo, bchi) bind(c,name='amrex_mlndlap_set_dot_mask')
    integer, dimension(1), intent(in) :: lo, hi, dlo, dhi, olo, ohi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(inout) :: dmsk(dlo(1):dhi(1))
    integer         , intent(in   ) :: omsk(olo(1):ohi(1))
  end subroutine amrex_mlndlap_set_dot_mask


  function amrex_mlndlap_any_fine_sync_cells (lo, hi, msk, mlo, mhi, fine_flag) result(r) &
       bind(c,name='amrex_mlndlap_any_fine_sync_cells')
    integer :: r
    integer, dimension(1), intent(in) :: lo, hi, mlo, mhi
    integer, intent(in   ) :: msk  ( mlo(1): mhi(1))
    integer, intent(in) :: fine_flag
  end function amrex_mlndlap_any_fine_sync_cells


  subroutine amrex_mlndlap_avgdown_coeff (lo, hi, crse, clo, chi, fine, flo, fhi, idim) &
       bind(c,name='amrex_mlndlap_avgdown_coeff')
    integer, dimension(1), intent(in) :: lo, hi, clo, chi, flo, fhi
    integer, intent(in) :: idim
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1))
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1))
  end subroutine amrex_mlndlap_avgdown_coeff


  subroutine amrex_mlndlap_fillbc_cc (sigma, slo, shi, dlo, dhi, bclo, bchi) &
       bind(c, name='amrex_mlndlap_fillbc_cc')
    integer, dimension(1), intent(in) :: slo, shi, dlo, dhi, bclo, bchi
    real(amrex_real), intent(inout) :: sigma(slo(1):shi(1))
  end subroutine amrex_mlndlap_fillbc_cc


  subroutine amrex_mlndlap_applybc (phi, hlo, hhi, dlo, dhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_applybc')
    integer, dimension(1) :: hlo, hhi, dlo, dhi, bclo, bchi
    real(amrex_real), intent(inout) :: phi(hlo(1):hhi(1))
  end subroutine amrex_mlndlap_applybc


  subroutine amrex_mlndlap_adotx_ha (lo, hi, y, ylo, yhi, x, xlo, xhi, &
       sx, sxlo, sxhi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_adotx_ha')
    integer, dimension(1), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, sxlo, sxhi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) ::  y( ylo(1): yhi(1))
    real(amrex_real), intent(in   ) ::  x( xlo(1): xhi(1))
    real(amrex_real), intent(in   ) :: sx(sxlo(1):sxhi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_adotx_ha


  subroutine amrex_mlndlap_adotx_aa (lo, hi, y, ylo, yhi, x, xlo, xhi, &
       sig, slo, shi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_adotx_aa')
    integer, dimension(1), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, slo, shi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) ::   y(ylo(1):yhi(1))
    real(amrex_real), intent(in   ) ::   x(xlo(1):xhi(1))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_adotx_aa


  subroutine amrex_mlndlap_jacobi_ha (lo, hi, sol, slo, shi, Ax, alo, ahi, rhs, rlo, rhi, &
       sx, sxlo, sxhi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_jacobi_ha')
    integer, dimension(1),intent(in) :: lo,hi,slo,shi,alo,ahi,rlo,rhi,sxlo,sxhi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1))
    real(amrex_real), intent(in   ) :: Ax ( alo(1): ahi(1))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1))
    real(amrex_real), intent(in   ) :: sx (sxlo(1):sxhi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_jacobi_ha


  subroutine amrex_mlndlap_jacobi_aa (lo, hi, sol, slo, shi, Ax, alo, ahi, rhs, rlo, rhi, &
       sig, sglo, sghi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_jacobi_aa')
    integer, dimension(1),intent(in) :: lo,hi,slo,shi,alo,ahi,rlo,rhi,sglo,sghi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1))
    real(amrex_real), intent(in   ) :: Ax ( alo(1): ahi(1))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1))
    real(amrex_real), intent(in   ) :: sig(sglo(1):sghi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_jacobi_aa


  subroutine amrex_mlndlap_gauss_seidel_ha (lo, hi, sol, slo, shi, rhs, rlo, rhi, &
       sx, sxlo, sxhi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_gauss_seidel_ha')
    integer, dimension(1),intent(in) :: lo,hi,slo,shi,rlo,rhi,sxlo,sxhi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1))
    real(amrex_real), intent(in   ) :: sx (sxlo(1):sxhi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_gauss_seidel_ha


  subroutine amrex_mlndlap_gauss_seidel_aa (lo, hi, sol, slo, shi, rhs, rlo, rhi, &
       sig, sglo, sghi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_gauss_seidel_aa')
    integer, dimension(1),intent(in) :: lo,hi,slo,shi,rlo,rhi,sglo,sghi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1))
    real(amrex_real), intent(in   ) :: sig(sglo(1):sghi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_gauss_seidel_aa


  subroutine amrex_mlndlap_restriction (lo, hi, crse, clo, chi, fine, flo, fhi, msk, mlo, mhi, &
       domlo, domhi, bclo, bchi) bind(c,name='amrex_mlndlap_restriction')
    integer, dimension(1), intent(in) :: lo, hi, clo, chi, flo, fhi, mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1))
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_restriction


  subroutine amrex_mlndlap_interpolation_ha (clo, chi, fine, fflo, ffhi, crse, cflo, cfhi, &
       sigx, sxlo, sxhi, msk, mlo, mhi, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_interpolation_ha')
    integer, dimension(1), intent(in) :: clo,chi,fflo,ffhi,cflo,cfhi,sxlo,sxhi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in   ) :: crse(cflo(1):cfhi(1))
    real(amrex_real), intent(inout) :: fine(fflo(1):ffhi(1))
    real(amrex_real), intent(in   ) :: sigx(sxlo(1):sxhi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_interpolation_ha


  subroutine amrex_mlndlap_interpolation_aa (clo, chi, fine, fflo, ffhi, crse, cflo, cfhi, &
       sig, sglo, sghi, msk, mlo, mhi, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_interpolation_aa')
    integer, dimension(1), intent(in) :: clo,chi,fflo,ffhi,cflo,cfhi,sglo,sghi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in   ) :: crse(cflo(1):cfhi(1))
    real(amrex_real), intent(inout) :: fine(fflo(1):ffhi(1))
    real(amrex_real), intent(in   ) :: sig (sglo(1):sghi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_interpolation_aa


  subroutine amrex_mlndlap_divu (lo, hi, rhs, rlo, rhi, vel, vlo, vhi, msk, mlo, mhi, &
       dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu')
    integer, dimension(1), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, mlo, mhi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_divu


  subroutine amrex_mlndlap_add_rhcc (lo, hi, rhs, rlo, rhi, rhcc, clo, chi, msk, mlo, mhi) &
       bind(c,name='amrex_mlndlap_add_rhcc')
    integer, dimension(1) :: lo, hi, rlo, rhi, clo, chi, mlo, mhi
    real(amrex_real), intent(inout) :: rhs (rlo(1):rhi(1))
    real(amrex_real), intent(in   ) :: rhcc(clo(1):chi(1))
    integer,          intent(in   ) :: msk (mlo(1):mhi(1))
  end subroutine amrex_mlndlap_add_rhcc


  subroutine amrex_mlndlap_mknewu (lo, hi, u, ulo, uhi, p, plo, phi, sig, slo, shi, dxinv) &
       bind(c,name='amrex_mlndlap_mknewu')
    integer, dimension(1), intent(in) :: lo, hi, ulo, uhi, plo, phi, slo, shi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) ::   u(ulo(1):uhi(1))
    real(amrex_real), intent(in   ) ::   p(plo(1):phi(1))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1))
  end subroutine amrex_mlndlap_mknewu


  subroutine amrex_mlndlap_divu_fine_contrib (clo, chi, cglo, cghi, rhs, rlo, rhi, &
       vel, vlo, vhi, frh, flo, fhi, msk, mlo, mhi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu_fine_contrib')
    integer, dimension(1), intent(in) :: clo, chi, cglo, cghi, rlo, rhi, vlo, vhi, &
         flo, fhi, mlo, mhi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1))
    real(amrex_real), intent(inout) :: frh(flo(1):fhi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_divu_fine_contrib


  subroutine amrex_mlndlap_divu_cf_contrib (lo, hi,  rhs, rlo, rhi, vel, vlo, vhi, dmsk, mlo, mhi, &
       ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi, fc, clo, chi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu_cf_contrib')
    integer, dimension(1), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, mlo, mhi, &
         nmlo, nmhi, cmlo, cmhi, clo, chi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1))
    real(amrex_real), intent(in   ) :: fc (clo(1):chi(1))
    integer, intent(in) :: dmsk(mlo(1):mhi(1))
    integer, intent(in) :: ndmsk(nmlo(1):nmhi(1))
    integer, intent(in) :: ccmsk(cmlo(1):cmhi(1))
  end subroutine amrex_mlndlap_divu_cf_contrib


  subroutine amrex_mlndlap_crse_resid (lo, hi, resid, rslo, rshi, rhs, rhlo, rhhi, msk, mlo, mhi) &
       bind(c, name='amrex_mlndlap_crse_resid')
    integer, dimension(1), intent(in) :: lo, hi, rslo, rshi, rhlo, rhhi, mlo, mhi
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
       fc, clo, chi, dxinv) &
       bind(c,name='amrex_mlndlap_res_cf_contrib')
    integer, dimension(1), intent(in) :: lo, hi, rlo, rhi, phlo, phhi, rhlo, rhhi, slo, shi, &
         mlo, mhi, nmlo, nmhi, cmlo, cmhi, clo, chi
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
