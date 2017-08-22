module eb_diffusion_module
  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

  public :: eb_diff_mol_3d

contains

  subroutine eb_diff_mol_3d (q, qd_lo, qd_hi, &
                     lo, hi, dx, &
                     flux1, fd1_lo, fd1_hi, &
                     flux2, fd2_lo, fd2_hi, &
                     flux3, fd3_lo, fd3_hi)

    use mempool_module, only : amrex_allocate, amrex_deallocate
    use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar, &
         qrho,qu,qv,qw,qp,qc,qeint,qtemp,qvar

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

  end subroutine eb_diff_mol_3d

end module eb_diffusion_module
