
subroutine amrex_eb_avgdown_sv (lo, hi, fine, flo, fhi, crse, clo, chi, flag, fglo, fghi, &
     fv, fvlo, fvhi, vfrc, vflo, vfhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown_sv')
  use amrex_fort_module, only : amrex_real
  use amrex_ebcellflag_module, only : is_covered_cell
  implicit none
  integer, intent(in) :: lo(3), hi(3), flo(3), fhi(3), clo(3), chi(3), fglo(3), fghi(3), &
       fvlo(3), fvhi(3), vflo(3), vfhi(3), lrat(3), ncomp
  real(amrex_real), intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),ncomp)
  real(amrex_real), intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3),ncomp)
  real(amrex_real), intent(in   ) :: fv  (fvlo(1):fvhi(1),fvlo(2):fvhi(2),fvlo(3):fvhi(3))
  real(amrex_real), intent(in   ) :: vfrc(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))
  integer         , intent(in   ) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

  integer :: i, j, k, ii, jj, kk, n, iref, jref, kref
  real(amrex_real) :: cv

  do n = 1, ncomp
     do k        = lo(3), hi(3)
        kk       = k * lrat(3)
        do j     = lo(2), hi(2)
           jj    = j * lrat(2)
           do i  = lo(1), hi(1)
              ii = i * lrat(1)
              crse(i,j,k,n) = 0.d0
              cv            = 0.d0
              if (is_covered_cell(flag(i,j,k))) then
                 do       kref = 0, lrat(3)-1
                    do    jref = 0, lrat(2)-1
                       do iref = 0, lrat(1)-1
                          cv = cv + fv(ii+iref,jj+jref,kk+kref)
                          crse(i,j,k,n) = crse(i,j,k,n) + &
                               fine(ii+iref,jj+jref,kk+kref,n)*fv(ii+iref,jj+jref,kk+kref)
                       end do
                    end do
                 end do
              else
                 do       kref = 0, lrat(3)-1
                    do    jref = 0, lrat(2)-1
                       do iref = 0, lrat(1)-1
                          cv = cv + (fv(ii+iref,jj+jref,kk+kref)*vfrc(ii+iref,jj+jref,kk+kref))
                          crse(i,j,k,n) = crse(i,j,k,n) + &
                               fine(ii+iref,jj+jref,kk+kref,n)*(fv(ii+iref,jj+jref,kk+kref)*vfrc(ii+iref,jj+jref,kk+kref))
                       end do
                    end do
                 end do
              endif
              crse(i,j,k,n) = crse(i,j,k,n) / cv
           end do
        end do
     end do
  end do

end subroutine amrex_eb_avgdown_sv
