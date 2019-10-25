module amrex_mlnodelap_1d_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use amrex_constants_module
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann, amrex_lo_inflow, amrex_lo_periodic
  implicit none

  ! external dirichlet at physical boundary or internal dirichlet at crse/fine boundary
  integer, parameter :: dirichlet = 1

  integer, parameter :: crse_cell = 0
  integer, parameter :: fine_cell = 1
  integer, parameter :: crse_node = 0
  integer, parameter :: crse_fine_node = 1
  integer, parameter :: fine_node = 2

  real(amrex_real), private, parameter :: eps = 1.d-100

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
       amrex_mlndlap_divu_fine_contrib, amrex_mlndlap_divu_cf_contrib, &
       amrex_mlndlap_rhcc_fine_contrib, amrex_mlndlap_rhcc_crse_contrib, &
       amrex_mlndlap_vel_cc_to_ct, amrex_mlndlap_mknewu_eb, &
       ! residual
       amrex_mlndlap_crse_resid, &
       amrex_mlndlap_res_fine_contrib, amrex_mlndlap_res_cf_contrib, &
       ! sync residual
       amrex_mlndlap_zero_fine

  ! RAP
  public:: &
       amrex_mlndlap_stencil_rap

#ifdef AMREX_USE_EB
  public:: amrex_mlndlap_set_integral, amrex_mlndlap_set_integral_eb, &
       amrex_mlndlap_set_connection, amrex_mlndlap_set_stencil_eb, &
       amrex_mlndlap_divu_eb, amrex_mlndlap_mknewu_eb
#endif

contains

  function amrex_mlndlap_any_fine_sync_cells (lo, hi, msk, mlo, mhi, fine_flag) result(r) &
       bind(c,name='amrex_mlndlap_any_fine_sync_cells')
    integer :: r
    integer, dimension(1), intent(in) :: lo, hi, mlo, mhi
    integer, intent(in   ) :: msk  ( mlo(1): mhi(1))
    integer, intent(in) :: fine_flag
  end function amrex_mlndlap_any_fine_sync_cells


  subroutine amrex_mlndlap_vel_cc_to_ct (lo, hi, vel, vlo, vhi, ovel, olo, ohi, vfrac, flo, fhi, &
       cent, clo, chi, flag, glo, ghi) bind(c,name='amrex_mlndlap_vel_cc_to_ct')
    integer, dimension(1), intent(in) :: lo, hi, vlo, vhi, olo, ohi, flo, fhi, clo, chi, glo, ghi
    real(amrex_real), intent(inout) ::   vel(vlo(1):vhi(1))
    real(amrex_real), intent(in   ) ::  ovel(olo(1):ohi(1))
    real(amrex_real), intent(in   ) :: vfrac(flo(1):fhi(1))
    real(amrex_real), intent(in   ) ::  cent(clo(1):chi(1))
    integer         , intent(in   ) ::  flag(glo(1):ghi(1))
  end subroutine amrex_mlndlap_vel_cc_to_ct


  subroutine amrex_mlndlap_mknewu_eb (lo, hi, u, ulo, uhi, p, plo, phi, sig, slo, shi, &
       vfrac, vlo, vhi, dxinv) bind(c,name='amrex_mlndlap_mknewu_eb')
    integer, dimension(1), intent(in) :: lo, hi, ulo, uhi, plo, phi, slo, shi, vlo, vhi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) ::   u(ulo(1):uhi(1))
    real(amrex_real), intent(in   ) ::   p(plo(1):phi(1))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1))
    real(amrex_real), intent(in   ) :: vfrac(vlo(1):vhi(1))
  end subroutine amrex_mlndlap_mknewu_eb


  subroutine amrex_mlndlap_divu_fine_contrib (clo, chi, cglo, cghi, rhs, rlo, rhi, &
       vel, vlo, vhi, frh, flo, fhi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_divu_fine_contrib')
    integer, dimension(1), intent(in) :: clo, chi, cglo, cghi, rlo, rhi, vlo, vhi, &
         flo, fhi, mlo, mhi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1))
    real(amrex_real), intent(inout) :: frh(flo(1):fhi(1))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1))
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


  subroutine amrex_mlndlap_rhcc_fine_contrib (clo, chi, cglo, cghi, rhs, rlo, rhi, &
       cc, cclo, cchi, msk, mlo, mhi) bind(c,name='amrex_mlndlap_rhcc_fine_contrib')
    integer, dimension(2), intent(in) :: clo, chi, cglo, cghi, rlo, rhi, cclo, cchi, mlo, mhi
    real(amrex_real), intent(inout) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(in   ) :: cc (cclo(1):cchi(1),cclo(2):cchi(2))
    integer         , intent(in   ) :: msk( mlo(1): mhi(1), mlo(2): mhi(2))
  end subroutine amrex_mlndlap_rhcc_fine_contrib


  subroutine amrex_mlndlap_rhcc_crse_contrib (lo, hi, crhs, rlo, rhi, rhcc, clo, chi, &
       dmsk, mlo, mhi, ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi) &
       bind(c,name='amrex_mlndlap_rhcc_crse_contrib')
    integer, dimension(2), intent(in) :: lo, hi, rlo, rhi, clo, chi, mlo, mhi, &
         nmlo, nmhi, cmlo, cmhi
    real(amrex_real), intent(inout) ::  crhs(rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(in   ) ::  rhcc(clo(1):chi(1),clo(2):chi(2))
    integer         , intent(in   ) ::  dmsk( mlo(1): mhi(1), mlo(2): mhi(2))
    integer         , intent(in   ) :: ndmsk(nmlo(1):nmhi(1),nmlo(2):nmhi(2))
    integer         , intent(in   ) :: ccmsk(cmlo(1):cmhi(1),cmlo(2):cmhi(2))
  end subroutine amrex_mlndlap_rhcc_crse_contrib


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


  subroutine amrex_mlndlap_stencil_rap (lo, hi, csten, clo, chi, fsten, flo, fhi) &
       bind(c,name='amrex_mlndlap_stencil_rap')
    integer, dimension(1), intent(in) :: lo, hi, clo, chi, flo, fhi
    real(amrex_real), intent(inout) :: csten(clo(1):chi(1),3)
    real(amrex_real), intent(in   ) :: fsten(flo(1):fhi(1),3)
  end subroutine amrex_mlndlap_stencil_rap


#ifdef AMREX_USE_EB

  subroutine amrex_mlndlap_set_integral (lo, hi, intg, glo, ghi) &
       bind(c,name='amrex_mlndlap_set_integral')
    integer, dimension(1) :: lo, hi, glo, ghi
    real(amrex_real), intent(inout) :: intg( glo(1): ghi(1),1)
  end subroutine amrex_mlndlap_set_integral

  subroutine amrex_mlndlap_set_integral_eb (lo, hi, intg, glo, ghi, flag, flo, fhi, &
       vol, vlo, vhi, ax, axlo, axhi, bcen, blo, bhi) &
       bind(c,name='amrex_mlndlap_set_integral_eb')
    use amrex_ebcellflag_module, only : is_single_valued_cell, is_regular_cell, is_covered_cell
    integer, dimension(1) :: lo, hi, glo, ghi, flo, fhi, axlo, vlo, vhi, axhi, blo, bhi
    real(amrex_real), intent(inout) :: intg( glo(1): ghi(1),1)
    real(amrex_real), intent(in   ) :: vol ( vlo(1): vhi(1))
    real(amrex_real), intent(in   ) :: ax  (axlo(1):axhi(1))
    real(amrex_real), intent(in   ) :: bcen( blo(1): bhi(1))
    integer         , intent(in   ) :: flag( flo(1): fhi(1))
  end subroutine amrex_mlndlap_set_integral_eb

  subroutine amrex_mlndlap_set_connection (lo, hi, conn, clo, chi, intg, glo, ghi, flag, flo, fhi, &
       vol, vlo, vhi) bind(c,name='amrex_mlndlap_set_connection')
    use amrex_ebcellflag_module, only : is_single_valued_cell, is_regular_cell, is_covered_cell
    integer, dimension(1) :: lo, hi, clo, chi, glo, ghi, flo, fhi, axlo, vlo, vhi
    real(amrex_real), intent(inout) :: conn( clo(1): chi(1),2)
    real(amrex_real), intent(inout) :: intg( glo(1): ghi(1),1)
    real(amrex_real), intent(in   ) :: vol ( vlo(1): vhi(1))
    integer         , intent(in   ) :: flag( flo(1): fhi(1))
  end subroutine amrex_mlndlap_set_connection


  subroutine amrex_mlndlap_set_stencil_eb (lo, hi, sten, tlo, thi, sigma, glo, ghi, &
       conn, clo, chi, dxinv) bind(c,name='amrex_mlndlap_set_stencil_eb')
    integer, dimension(1), intent(in) :: lo, hi, tlo, thi, glo, ghi, clo, chi
    real(amrex_real), intent(inout) ::  sten(tlo(1):thi(1),3)
    real(amrex_real), intent(in   ) :: sigma(glo(1):ghi(1))
    real(amrex_real), intent(in   ) ::  conn(clo(1):chi(1),2)
    real(amrex_real), intent(in) :: dxinv(1)
  end subroutine amrex_mlndlap_set_stencil_eb


  subroutine amrex_mlndlap_divu_eb (lo, hi, rhs, rlo, rhi, vel, vlo, vhi, vfrac, flo, fhi, &
       intg, glo, ghi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_divu_eb')
    integer, dimension(1), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, flo, fhi, glo, ghi, mlo, mhi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1))
    real(amrex_real), intent(in   ) :: vfrac(flo(1):fhi(1))
    real(amrex_real), intent(in   ) :: intg(glo(1):ghi(1))
    integer, intent(in) :: msk(mlo(1):mhi(1))
  end subroutine amrex_mlndlap_divu_eb


  subroutine amrex_mlndlap_mknewu_eb (lo, hi, u, ulo, uhi, p, plo, phi, sig, slo, shi, &
       vfrac, vlo, vhi, intg, glo, ghi, dxinv) bind(c,name='amrex_mlndlap_mknewu_eb')
    integer, dimension(1), intent(in) :: lo, hi, ulo, uhi, plo, phi, slo, shi, vlo, vhi, glo, ghi
    real(amrex_real), intent(in) :: dxinv(1)
    real(amrex_real), intent(inout) ::   u(ulo(1):uhi(1))
    real(amrex_real), intent(in   ) ::   p(plo(1):phi(1))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1))
    real(amrex_real), intent(in   )::vfrac(vlo(1):vhi(1))
    real(amrex_real), intent(in   ) ::intg(glo(1):ghi(1))
  end subroutine amrex_mlndlap_mknewu_eb

#endif

end module amrex_mlnodelap_1d_module
