
subroutine amrex_eb_avgdown_sv (lo, hi, fine, flo, fhi, crse, clo, chi, flag, fglo, fghi, &
     fv, fvlo, fvhi, vfrc, vflo, vfhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown_sv')
  use amrex_fort_module, only : amrex_real
  use amrex_ebcellflag_module, only : is_covered_cell
  implicit none
  integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), clo(2), chi(2), fglo(2), fghi(2), &
       fvlo(2), fvhi(2), vflo(2), vfhi(2), lrat(2), ncomp
  real(amrex_real), intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2),ncomp)
  real(amrex_real), intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2),ncomp)
  real(amrex_real), intent(in   ) :: fv  (fvlo(1):fvhi(1),fvlo(2):fvhi(2))
  real(amrex_real), intent(in   ) :: vfrc(vflo(1):vfhi(1),vflo(2):vfhi(2))
  integer         , intent(in   ) :: flag(fglo(1):fghi(1),fglo(2):fghi(2))

  integer :: i, j, ii, jj, n, iref, jref
  real(amrex_real) :: cv

  do n = 1, ncomp
     do j     = lo(2), hi(2)
        jj    = j * lrat(2)
        do i  = lo(1), hi(1)
           ii = i * lrat(1)
           crse(i,j,n) = 0.d0
           cv          = 0.d0
           if (is_covered_cell(flag(i,j))) then
              do    jref = 0, lrat(2)-1
                 do iref = 0, lrat(1)-1
                    cv          = cv          +                         fv(ii+iref,jj+jref)
                    crse(i,j,n) = crse(i,j,n) + fine(ii+iref,jj+jref,n)*fv(ii+iref,jj+jref)
                 end do
              end do
           else
              do    jref = 0, lrat(2)-1
                 do iref = 0, lrat(1)-1
                    cv          = cv          +                         (fv(ii+iref,jj+jref) &
                          * vfrc(ii+iref,jj+jref))
                    crse(i,j,n) = crse(i,j,n) + fine(ii+iref,jj+jref,n)*(fv(ii+iref,jj+jref) &
                         * vfrc(ii+iref,jj+jref))
                 end do
              end do
           endif
           crse(i,j,n) = crse(i,j,n) / cv
        end do
     end do
  end do

end subroutine amrex_eb_avgdown_sv
