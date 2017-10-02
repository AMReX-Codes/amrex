module amrex_ebinterp_module

  use amrex_fort_module, only : amrex_real
  implicit none
  private
  public :: amrex_ebinterp_pc_sv

contains

  subroutine amrex_ebinterp_pc_sv (tflo, tfhi, tclo, tchi, crse, clo, chi, fine, flo, fhi, &
       ncomp, ratio, cdomainlo, cdomainhi, cflag, cflo, cfhi) &
       bind(c,name='amrex_ebinterp_pc_sv')
    integer, intent(in) :: tclo(1), tchi(1), tflo(1), tfhi(1), clo(1), chi(1), flo(1), fhi(1)
    integer, intent(in) :: ncomp, ratio(1), cflo(1), cfhi(1)
    integer, intent(in) :: cdomainlo(1), cdomainhi(1)
    real(amrex_real), intent(in   ) :: crse(clo(1):chi(1),ncomp)
    real(amrex_real), intent(inout) :: fine(flo(1):fhi(1),ncomp)
    integer, intent(in) :: cflag(cflo(1):cfhi(1))

    stop '1d amrex_ebinterp_pc_sv not support'

  end subroutine amrex_ebinterp_pc_sv

end module amrex_ebinterp_module
