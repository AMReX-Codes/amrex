module amrex_ebinterp_module

  use amrex_fort_module, only : amrex_real
  use amrex_ebcellflag_module, only : num_neighbor_cells
  implicit none
  private
  public :: amrex_ebinterp_pc_sv

contains

  subroutine amrex_ebinterp_pc_sv (tflo, tfhi, tclo, tchi, crse, clo, chi, fine, flo, fhi, &
       ncomp, ratio, cdomainlo, cdomainhi, cflag, cflo, cfhi) &
       bind(c,name='amrex_ebinterp_pc_sv')
    integer, intent(in) :: tclo(2), tchi(2), tflo(2), tfhi(2), clo(2), chi(2), flo(2), fhi(2)
    integer, intent(in) :: ncomp, ratio(2), cflo(2), cfhi(2)
    integer, intent(in) :: cdomainlo(2), cdomainhi(2)
    real(amrex_real), intent(in   ) :: crse(clo(1):chi(1),clo(2):chi(2),ncomp)
    real(amrex_real), intent(inout) :: fine(flo(1):fhi(1),flo(2):fhi(2),ncomp)
    integer, intent(in) :: cflag(cflo(1):cfhi(1),cflo(2):cfhi(2))

    integer :: ic, jc, n, i, j, jmin, jmax, imin, imax, glo(2), ghi(2)
    integer :: num_expected_ngbrs

    integer :: ngbr(-1:1,-1:1)

    do jc = tclo(2)+1, tchi(2)-1
       jmin = max(jc*ratio(2)           ,tflo(2))
       jmax = min(jc*ratio(2)+ratio(2)-1,tfhi(2))

       glo(2) = max(jc-1, cdomainlo(2))
       ghi(2) = min(jc+1, cdomainhi(2))

       do ic = tclo(1)+1, tchi(1)-1
          imin = max(ic*ratio(1)           ,tflo(1))
          imax = min(ic*ratio(1)+ratio(1)-1,tfhi(1))

          glo(1) = max(ic-1, cdomainlo(1))
          ghi(1) = min(ic+1, cdomainhi(1))

          num_expected_ngbrs = (ghi(1)-glo(1)+1)*(ghi(2)-glo(2)+1)

          if (num_neighbor_cells(cflag(ic,jc)) .lt. num_expected_ngbrs) then

             do n = 1, ncomp
                do    j = jmin, jmax
                   do i = imin, imax
                      fine(i,j,n) = crse(ic,jc,n)
                   end do
                end do
             end do
          end if

       end do
    end do

  end subroutine amrex_ebinterp_pc_sv

end module amrex_ebinterp_module
