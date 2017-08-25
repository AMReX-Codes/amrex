module diffusion_module
  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

  public :: diff_mol_3d

contains

  subroutine diff_mol_3d (q, qd_lo, qd_hi, &
                     lo, hi, dx, &
                     flux1, fd1_lo, fd1_hi, &
                     flux2, fd2_lo, fd2_hi, &
                     flux3, fd3_lo, fd3_hi)

    use mempool_module, only : amrex_allocate, amrex_deallocate
    use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar, &
         qrho,qu,qv,qw,qp,qc,qeint,qtemp,qvar
    use diff_coef_module, only : compute_diff_coef

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: fd1_lo(3), fd1_hi(3)
    integer, intent(in) :: fd2_lo(3), fd2_hi(3)
    integer, intent(in) :: fd3_lo(3), fd3_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in   ) ::     q( qd_lo(1): qd_hi(1), qd_lo(2): qd_hi(2), qd_lo(3): qd_hi(3),QVAR)
    real(rt), intent(inout) :: flux1(fd1_lo(1):fd1_hi(1),fd1_lo(2):fd1_hi(2),fd1_lo(3):fd1_hi(3),NVAR)
    real(rt), intent(inout) :: flux2(fd2_lo(1):fd2_hi(1),fd2_lo(2):fd2_hi(2),fd2_lo(3):fd2_hi(3),NVAR)
    real(rt), intent(inout) :: flux3(fd3_lo(1):fd3_hi(1),fd3_lo(2):fd3_hi(2),fd3_lo(3):fd3_hi(3),NVAR)

    integer :: i,j,k, clo(3), chi(3)
    integer, parameter :: nextra = 0
    real(rt), dimension(:,:,:), pointer, contiguous :: lambda, mu, xi
    
    clo = lo-nextra-1
    chi = lo+nextra+1

    call amrex_allocate(lambda, clo, chi)
    call amrex_allocate(mu, clo, chi)
    call amrex_allocate(xi, clo, chi)

    call compute_diff_coef(q, qd_lo, qd_hi, lambda, mu, xi, clo, chi)

    call amrex_deallocate(lambda)
    call amrex_deallocate(mu)
    call amrex_deallocate(xi)

  end subroutine diff_mol_3d

end module diffusion_module
