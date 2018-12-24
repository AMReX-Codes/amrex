module amrex_eb_bc_fill_module

! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
! for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains

  subroutine amrex_eb_phifill(phi,phi_lo,phi_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="amrex_eb_phifill")

    use amrex_fort_module, only : bl_spacedim, amrex_real
    use amrex_filcc_module, only : amrex_filccn

    implicit none

    integer      :: phi_lo(3),phi_hi(3)
    integer      :: bc(bl_spacedim,2)
    integer      :: domlo(3), domhi(3)
    real(amrex_real) :: delta(3), xlo(3), time
    real(amrex_real) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

    call amrex_filccn(phi_lo, phi_hi, phi, phi_lo, phi_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine amrex_eb_phifill

end module amrex_eb_bc_fill_module
