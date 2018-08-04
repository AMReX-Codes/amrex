module amrex_eb_avgdown_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_eb_avgdown_sv, amrex_eb_avgdown, amrex_eb_avgdown_faces

contains

  subroutine amrex_eb_avgdown_sv (lo, hi, fine, flo, fhi, crse, clo, chi, &
       fv, fvlo, fvhi, vfrc, vflo, vfhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown_sv')
    implicit none
    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), clo(2), chi(2), &
         fvlo(2), fvhi(2), vflo(2), vfhi(2), lrat(2), ncomp
    real(amrex_real), intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2),ncomp)
    real(amrex_real), intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2),ncomp)
    real(amrex_real), intent(in   ) :: fv  (fvlo(1):fvhi(1),fvlo(2):fvhi(2))
    real(amrex_real), intent(in   ) :: vfrc(vflo(1):vfhi(1),vflo(2):vfhi(2))
    
    integer :: i, j, ii, jj, n, iref, jref
    real(amrex_real) :: cv
    
    do n = 1, ncomp
       do j     = lo(2), hi(2)
          jj    = j * lrat(2)
          do i  = lo(1), hi(1)
             ii = i * lrat(1)
             crse(i,j,n) = 0.d0
             cv          = 0.d0
             do    jref = 0, lrat(2)-1
                do iref = 0, lrat(1)-1
                   cv          = cv          +                         (fv(ii+iref,jj+jref) &
                        * vfrc(ii+iref,jj+jref))
                   crse(i,j,n) = crse(i,j,n) + fine(ii+iref,jj+jref,n)*(fv(ii+iref,jj+jref) &
                        * vfrc(ii+iref,jj+jref))
                end do
             end do
             if (cv .gt. 1.d-30) then
                crse(i,j,n) = crse(i,j,n) / cv
             else
                crse(i,j,n) = fine(ii,jj,n)  ! covered cell
             end if
          end do
       end do
    end do
  end subroutine amrex_eb_avgdown_sv


  subroutine amrex_eb_avgdown (lo, hi, fine, flo, fhi, crse, clo, chi, &
       vfrc, vflo, vfhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown')
    implicit none
    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), clo(2), chi(2), &
         vflo(2), vfhi(2), lrat(2), ncomp
    real(amrex_real), intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2),ncomp)
    real(amrex_real), intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2),ncomp)
    real(amrex_real), intent(in   ) :: vfrc(vflo(1):vfhi(1),vflo(2):vfhi(2))
    
    integer :: i, j, ii, jj, n, iref, jref
    real(amrex_real) :: cv
    
    do n = 1, ncomp
       do j     = lo(2), hi(2)
          jj    = j * lrat(2)
          do i  = lo(1), hi(1)
             ii = i * lrat(1)
             crse(i,j,n) = 0.d0
             cv          = 0.d0
             do    jref = 0, lrat(2)-1
                do iref = 0, lrat(1)-1
                   cv  = cv + vfrc(ii+iref,jj+jref)
                   crse(i,j,n) = crse(i,j,n) + fine(ii+iref,jj+jref,n) * vfrc(ii+iref,jj+jref)
                end do
             end do
             if (cv .gt. 1.d-30) then
                crse(i,j,n) = crse(i,j,n) / cv
             else
                crse(i,j,n) = fine(ii,jj,n)  ! covered cell
             end if
          end do
       end do
    end do
  end subroutine amrex_eb_avgdown

  subroutine amrex_eb_avgdown_faces (lo, hi, fine, flo, fhi, crse, clo, chi, &
       ap, aplo, aphi, lrat, idir, ncomp) bind(c,name='amrex_eb_avgdown_faces')
    implicit none
    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), clo(2), chi(2), aplo(2), aphi(2), lrat(2), idir ,ncomp 
    real(amrex_real), intent(in   )   :: fine( flo(1): fhi(1), flo(2): fhi(2), ncomp)
    real(amrex_real), intent(inout)   :: crse( clo(1): chi(1), clo(2): chi(2), ncomp) 
    real(amrex_real), intent(in   )   ::  ap(aplo(1):aphi(1),aplo(2):aphi(2)) 

    integer :: i, j, ii, jj, n, iref, jref, facx, facy
    real(amrex_real) :: fa 
  
    facx = lrat(1) 
    facy = lrat(2)
    if(idir == 0) then
      do n = 1, ncomp 
         do j       = lo(2), hi(2) 
            jj      = j*facy
            do i    = lo(1), hi(1)
               ii   = i*facx
               crse(i,j,n) = 0.d0 
               fa          = 0.d0 
               do    jref  = 0, facy-1
                     fa          = fa + ap(ii, jj+jref) 
                     crse(i,j,n) = crse(i,j,n) + ap(ii,jj+jref)*fine(ii,jj+jref,n)
               enddo
               if(fa.gt.1.d-30) then 
                 crse(i,j,n) = crse(i,j,n)/fa
               else 
                 crse(i,j,n) = fine(ii,jj,n) !covered face
               endif 
            enddo
         enddo 
      enddo
    else
      do n = 1, ncomp
         do j       = lo(2), hi(2)
            jj      = j*facy
            do i    = lo(1), hi(1) 
               ii   = i*facx
               crse(i,j,n) = 0.d0
               fa          = 0.d0 
               do    iref  = 0, facx-1
                     fa          = fa + ap(ii+iref, jj) 
                     crse(i,j,n) = crse(i,j,n) + ap(ii+iref,jj)*fine(ii+iref,jj,n)
               enddo
               if(fa.gt.1.d-30) then 
                 crse(i,j,n) = crse(i,j,n)/fa
               else
                 crse(i,j,n) = fine(ii,jj,n) !covered face 
               endif
             enddo
         enddo
      enddo 
   endif
  end subroutine amrex_eb_avgdown_faces 

end module amrex_eb_avgdown_module
