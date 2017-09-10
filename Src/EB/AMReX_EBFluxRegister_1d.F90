module amrex_eb_flux_reg_1d_module

  implicit none
  private

  public :: amrex_eb_flux_reg_crseadd

contains

  subroutine amrex_eb_flux_reg_crseadd (lo, hi, fx, d, dlo, dhi, fxlo, fxhi, dx, dt, nc) &
       bind(c,name='amrex_eb_flux_reg_crseadd')
    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    integer :: lo, hi, dlo, dhi, fxlo, fxhi, nc
    real(rt), intent(in) :: d(dlo:dhi,nc)
    real(rt), intent(in) :: fx(fxlo:fxhi,nc)
    real(rt), intent(in) :: dx, dt
    
    call amrex_error("EB not supported for 1D")
  end subroutine amrex_eb_flux_reg_crseadd

end module amrex_eb_flux_reg_1d_module
