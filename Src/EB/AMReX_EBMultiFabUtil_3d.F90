module amrex_eb_avgdown_module
  use amrex_fort_module, only : amrex_real
  implicit none
  
  private
  public :: amrex_eb_avgdown_sv, amrex_eb_avgdown, amrex_eb_avgdown_faces

contains

  subroutine amrex_eb_avgdown_sv (lo, hi, fine, flo, fhi, crse, clo, chi, &
       fv, fvlo, fvhi, vfrc, vflo, vfhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown_sv')
    integer, intent(in) :: lo(3), hi(3), flo(3), fhi(3), clo(3), chi(3), &
         fvlo(3), fvhi(3), vflo(3), vfhi(3), lrat(3), ncomp
    real(amrex_real), intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),ncomp)
    real(amrex_real), intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3),ncomp)
    real(amrex_real), intent(in   ) :: fv  (fvlo(1):fvhi(1),fvlo(2):fvhi(2),fvlo(3):fvhi(3))
    real(amrex_real), intent(in   ) :: vfrc(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))
    
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
                do       kref = 0, lrat(3)-1
                   do    jref = 0, lrat(2)-1
                      do iref = 0, lrat(1)-1
                         cv = cv + (fv(ii+iref,jj+jref,kk+kref)*vfrc(ii+iref,jj+jref,kk+kref))
                         crse(i,j,k,n) = crse(i,j,k,n) + &
                              fine(ii+iref,jj+jref,kk+kref,n)*(fv(ii+iref,jj+jref,kk+kref)*vfrc(ii+iref,jj+jref,kk+kref))
                      end do
                   end do
                end do
                if (cv .gt. 1.d-30) then
                   crse(i,j,k,n) = crse(i,j,k,n) / cv
                else
                   crse(i,j,k,n) = fine(ii,jj,kk,n)  ! covered cell
                end if
             end do
          end do
       end do
    end do
  end subroutine amrex_eb_avgdown_sv


  subroutine amrex_eb_avgdown (lo, hi, fine, flo, fhi, crse, clo, chi, &
       vfrc, vflo, vfhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown')
    integer, intent(in) :: lo(3), hi(3), flo(3), fhi(3), clo(3), chi(3), &
         vflo(3), vfhi(3), lrat(3), ncomp
    real(amrex_real), intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),ncomp)
    real(amrex_real), intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3),ncomp)
    real(amrex_real), intent(in   ) :: vfrc(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))
    
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
                do       kref = 0, lrat(3)-1
                   do    jref = 0, lrat(2)-1
                      do iref = 0, lrat(1)-1
                         cv = cv + vfrc(ii+iref,jj+jref,kk+kref)
                         crse(i,j,k,n) = crse(i,j,k,n) + &
                              fine(ii+iref,jj+jref,kk+kref,n)*vfrc(ii+iref,jj+jref,kk+kref)
                      end do
                   end do
                end do
                if (cv .gt. 1.d-30) then
                   crse(i,j,k,n) = crse(i,j,k,n) / cv
                else
                   crse(i,j,k,n) = fine(ii,jj,kk,n)  ! covered cell
                end if
             end do
          end do
       end do
    end do
  end subroutine amrex_eb_avgdown
  
  subroutine amrex_eb_avgdown_faces (lo, hi, fine, flo, fhi, crse, clo, chi, &
       ap, aplo, aphi, lrat, idir, ncomp) bind(c,name='amrex_eb_avgdown_faces')
    integer, dimension(3), intent(in) :: lo, hi, flo, fhi, clo,chi, aplo, aphi, lrat
    integer,               intent(in) :: idir, ncomp
    real(amrex_real),   intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),ncomp) 
    real(amrex_real),   intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3),ncomp)
    real(amrex_real),   intent(in   ) ::   ap(aplo(1):aphi(1),aplo(2):aphi(2),aplo(3):aphi(3))

    integer  :: i, j, k, ii, jj, kk, n, iref, jref, kref
    real(amrex_real) :: fa 
 
    if(idir.eq.0) then 
      do n              = 1, ncomp
         do k           = lo(3), hi(3)
            kk          = k*lrat(3)
            do j        = lo(2), hi(2)
               jj       = j*lrat(2)
               do i     = lo(1), hi(1)
                  ii    = i*lrat(1)
                  crse(i,j,k,n) = 0.d0
                  fa            = 0.d0
                  do    kref    = 0, lrat(3)-1
                    do  jref    = 0, lrat(2)-1
                        fa            = fa + ap(ii,jj+jref,kk+kref)
                        crse(i,j,k,n) = crse(i,j,k,n) + ap(ii,jj+jref,kk+kref)*fine(ii,jj+jref,kk+kref,n)
                    enddo
                  enddo
                  if(fa.gt.1.d-30) then 
                    crse(i,j,k,n) = crse(i,j,k,n)/fa
                  else
                    crse(i,j,k,n) = fine(ii,jj,kk,n) !covered face
                  endif
               enddo
            enddo
         enddo
      enddo 
    elseif(idir.eq.1) then 
      do n             = 1, ncomp   
         do k          = lo(3), hi(3)
            kk         = k*lrat(3)
            do j       = lo(2), hi(2)
               jj      = j*lrat(2)
               do i    = lo(1), hi(1)
                  ii   = i*lrat(1)
                  crse(i,j,k,n) = 0.d0
                  fa            = 0.d0
                  do    kref    = 0, lrat(3)-1
                    do  iref    = 0, lrat(1)-1
                        fa            = fa + ap(ii+iref, jj, k+kref)
                        crse(i,j,k,n) = crse(i,j,k,n) + ap(ii+iref,jj,kk+kref)*fine(ii+iref,jj,kk+kref,n)
                    enddo
                  enddo
                  if(fa.gt.1.d-30) then
                    crse(i,j,k,n) = crse(i,j,k,n)/fa
                  else
                    crse(i,j,k,n) = fine(ii,jj,kk,n) !covered face
                  endif
               enddo
            enddo
         enddo
      enddo
    else
      do n            = 1, ncomp
         do k         = lo(3), hi(3)
            kk        = k*lrat(3)
            do j      = lo(2), hi(2)
               jj     = j*lrat(2)
               do i   = lo(1), hi(1)
                  ii  = i*lrat(1)
                  crse(i,j,k,n) = 0.d0
                  fa            = 0.d0
                  do    jref    = 0, lrat(2)-1
                    do  iref    = 0, lrat(1)-1
                        fa            = fa + ap(ii+iref,jj+jref,kk)
                        crse(i,j,k,n) = crse(i,j,k,n) + ap(ii+iref,jj+jref,kk)*fine(ii+iref,jj+jref,kk,n)
                    enddo
                  enddo
                  if(fa.gt.1.d0-30) then
                    crse(i,j,k,n) = crse(i,j,k,n)/fa
                  else
                    crse(i,j,k,n) = fine(ii,jj,kk,n) !covered face
                  endif
               enddo
            enddo
         enddo
      enddo
    endif
  end subroutine amrex_eb_avgdown_faces

end module amrex_eb_avgdown_module
