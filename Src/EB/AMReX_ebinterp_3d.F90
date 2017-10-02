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
    integer, intent(in) :: tclo(3), tchi(3), tflo(3), tfhi(3), clo(3), chi(3), flo(3), fhi(3)
    integer, intent(in) :: ncomp, ratio(3), cflo(3), cfhi(3)
    integer, intent(in) :: cdomainlo(3), cdomainhi(3)
    real(amrex_real), intent(in   ) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),ncomp)
    real(amrex_real), intent(inout) :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),ncomp)
    integer, intent(in) :: cflag(cflo(1):cfhi(1),cflo(2):cfhi(2),cflo(3):cfhi(3))

    integer :: ic, jc, kc, n, i, j, k, imin, imax, jmin, jmax, kmin, kmax, glo(3), ghi(3)
    integer :: num_expected_ngbrs

    do kc = tclo(3)+1, tchi(3)-1
       kmin = max(kc*ratio(3)           ,tflo(3))
       kmax = min(kc*ratio(3)+ratio(3)-1,tfhi(3))

       glo(3) = max(kc-1, cdomainlo(3))
       ghi(3) = min(kc+1, cdomainhi(3))

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
             
             num_expected_ngbrs = (ghi(1)-glo(1)+1)*(ghi(2)-glo(2)+1)*(ghi(3)-glo(3)+1)

             if (num_neighbor_cells(cflag(ic,jc,kc)) .lt. num_expected_ngbrs) then
                do n = 1, ncomp
                   do      k = kmin, kmax
                      do    j = jmin, jmax
                         do i = imin, imax
                            fine(i,j,k,n) = crse(ic,jc,kc,n)
                         end do
                      end do
                   end do
                end do
             end if
             
          end do
       end do
    end do

  end subroutine amrex_ebinterp_pc_sv

end module amrex_ebinterp_module
