module amrex_eb_flux_reg_3d_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_eb_flux_reg_nd_module, only : crse_cell, crse_fine_boundary_cell, fine_cell, &
       reredistribution_threshold
  implicit none
  private

  public :: amrex_eb_flux_reg_crseadd_va, amrex_eb_flux_reg_fineadd_va, &
       amrex_eb_flux_reg_fineadd_dm

contains

  subroutine amrex_eb_flux_reg_crseadd_va (lo, hi, d, dlo, dhi, flag, fglo, fghi, &
       fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, &
       vfrac, vlo, vhi, ax, axlo, axhi, ay, aylo, ayhi, az, azlo, azhi, &
       dx, dt, nc) &
       bind(c,name='amrex_eb_flux_reg_crseadd_va')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, fglo, fghi, fxlo, fxhi, fylo, fyhi, fzlo, fzhi,&
         vlo, vhi, axlo, axhi, aylo, ayhi, azlo, azhi
    integer, intent(in) :: nc
    integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    real(rt), intent(in   ) :: fx   (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),nc)
    real(rt), intent(in   ) :: fy   (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),nc)
    real(rt), intent(in   ) :: fz   (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),nc)
    real(rt), intent(in   ) :: vfrac( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3))
    real(rt), intent(in   ) :: ax   (axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(rt), intent(in   ) :: ay   (aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(rt), intent(in   ) :: az   (azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(rt), intent(inout) :: d ( dlo(1): dhi(1), dlo(2): dhi(2), dlo(3): dhi(3),nc)
    real(rt), intent(in) :: dx(3), dt

    integer :: i,j,k
    real(rt) :: dtdx, dtdy, dtdz, volinv

    dtdx = dt/dx(1)
    dtdy = dt/dx(2)
    dtdz = dt/dx(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (flag(i,j,k).eq.crse_fine_boundary_cell .and. vfrac(i,j,k).gt.1.d-14) then

                volinv = 1._rt/vfrac(i,j,k)

                if (flag(i-1,j,k).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) - dtdx*fx(i,j,k,:)*(ax(i,j,k)*volinv)
                else if (flag(i+1,j,k).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) + dtdx*fx(i+1,j,k,:)*(ax(i+1,j,k)*volinv)
                end if

                if (flag(i,j-1,k).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) - dtdy*fy(i,j,k,:)*(ay(i,j,k)*volinv)
                else if (flag(i,j+1,k).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) + dtdy*fy(i,j+1,k,:)*(ay(i,j+1,k)*volinv)
                end if

                if (flag(i,j,k-1).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) - dtdz*fz(i,j,k,:)*(az(i,j,k)*volinv)
                else if (flag(i,j,k+1).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) + dtdz*fz(i,j,k+1,:)*(az(i,j,k+1)*volinv)
                end if

             end if
          end do
       end do
    end do

  end subroutine amrex_eb_flux_reg_crseadd_va


  subroutine amrex_eb_flux_reg_fineadd_va (lo, hi, d, dlo, dhi, f, flo, fhi, &
       cvol, clo, chi, vfrac, vlo, vhi, ax, axlo, axhi, ay, aylo, ayhi, az, azlo, azhi, &
       dx, dt, nc, dir, side, ratio) &
       bind(c,name='amrex_eb_flux_reg_fineadd_va')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, flo, fhi, ratio, &
         vlo, vhi, clo, chi, axlo, axhi, aylo, ayhi, azlo, azhi
    integer, intent(in) :: nc, dir, side
    real(rt), intent(in) :: dx(3), dt
    real(rt), intent(in   ) :: f(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),nc)
    real(rt), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nc)
    real(rt)                :: cvol ( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    real(rt), intent(in   ) :: vfrac( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3))
    real(rt), intent(in   ) :: ax   (axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(rt), intent(in   ) :: ay   (aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(rt), intent(in   ) :: az   (azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))

    integer :: i,j,k,n,ii,jj,kk,ioff,joff,koff
    real(rt) :: fac

    ! dx is fine.  vfrac, ax, ay, and az are fine as well.
    ! lo and hi are also relative to the fine box.

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             cvol(i,j,k) = sum(vfrac(i*ratio(1):i*ratio(1)+ratio(1)-1,  &
                  &                  j*ratio(2):j*ratio(2)+ratio(2)-1, &
                  &                  k*ratio(3):k*ratio(3)+ratio(3)-1))
             if (cvol(i,j,k).gt.1.d-14) then
                cvol(i,j,k) = 1._rt/cvol(i,j,k)
             else
                cvol(i,j,k) = 0._rt
             end if
          end do
       end do
    end do

    if (dir .eq. 0) then ! x-direction, lo(1) == hi(1)

       fac = dt / dx(1)

       if (side .eq. 0) then ! lo-side

          i = lo(1)
          ii = (i+1)*ratio(1) !!!
          
          do n = 1, nc
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do koff = 0, ratio(3)-1
                      kk =  k*ratio(3)+koff
                      do joff = 0, ratio(2)-1
                         jj = j*ratio(2)+joff
                         d(i,j,k,n) = d(i,j,k,n) - fac*f(ii,jj,kk,n)*(ax(ii,jj,kk)*cvol(i,j,k))
                      end do
                   end do
                end do
             end do
          end do

       else ! hi-side

          i = lo(1)
          ii = i*ratio(1) !!!

          do n = 1, nc
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do koff = 0, ratio(3)-1
                      kk = k*ratio(3)+koff
                      do joff = 0, ratio(2)-1
                         jj = j*ratio(2)+joff
                         d(i,j,k,n) = d(i,j,k,n) + fac*f(ii,jj,kk,n)*(ax(ii,jj,kk)*cvol(i,j,k))
                      end do
                   end do
                end do
             end do
          end do

       end if

    else if (dir .eq. 1) then ! y-direction, lo(2) == hi(2)

       fac = dt / dx(2)

       if (side .eq. 0) then ! lo-side

          j = lo(2)
          jj = (j+1)*ratio(2) !!!

          do n = 1, nc
             do k = lo(3), hi(3)                
                do koff = 0, ratio(3)-1
                   kk = k*ratio(3)+koff
                   do ioff = 0, ratio(1)-1
                      do i = lo(1), hi(1)
                         ii = i*ratio(1)+ioff
                         d(i,j,k,n) = d(i,j,k,n) - fac*f(ii,jj,kk,n)*(ay(ii,jj,kk)*cvol(i,j,k))
                      end do
                   end do
                end do
             end do
          end do
             
       else ! hi-side

          j = lo(2)
          jj = j*ratio(2) !!!

          do n = 1, nc
             do k = lo(3), hi(3)
                do koff = 0, ratio(3)-1
                   kk = k*ratio(3)+koff
                   do ioff = 0, ratio(1)-1
                      do i = lo(1), hi(1)
                         ii = i*ratio(1)+ioff
                         d(i,j,k,n) = d(i,j,k,n) + fac*f(ii,jj,kk,n)*(ay(ii,jj,kk)*cvol(i,j,k))
                      end do
                   end do
                end do
             end do
          end do

       end if

    else  ! z-direction, lo(3) == hi(3)

       fac = dt / dx(3)

       if (side .eq. 0) then ! lo-side

          k = lo(3)
          kk = (k+1)*ratio(3) !!!

          do n = 1, nc
             do j = lo(2), hi(2)
                do joff = 0, ratio(2)-1
                   jj = j*ratio(2)+joff
                   do ioff = 0, ratio(1)-1
                      do i = lo(1), hi(1)
                         ii = i*ratio(1)+ioff
                         d(i,j,k,n) = d(i,j,k,n) - fac*f(ii,jj,kk,n)*(az(ii,jj,kk)*cvol(i,j,k))
                      end do
                   end do
                end do
             end do
          end do

       else ! hi-side

          k = lo(3)
          kk = k*ratio(3) !!!

          do n = 1, nc
             do j = lo(2), hi(2)
                do joff = 0, ratio(2)-1
                   jj = j*ratio(2)+joff
                   do ioff = 0, ratio(1)-1
                      do i = lo(1), hi(1)
                         ii = i*ratio(1)+ioff
                         d(i,j,k,n) = d(i,j,k,n) + fac*f(ii,jj,kk,n)*(az(ii,jj,kk)*cvol(i,j,k))
                      end do
                   end do
                end do
             end do
          end do

       end if
    end if
  end subroutine amrex_eb_flux_reg_fineadd_va


  subroutine amrex_eb_flux_reg_fineadd_dm (lo, hi, d, dlo, dhi, dm, mlo, mhi, &
       cvol, clo, chi, vfrac, vlo, vhi, dx, nc, ratio) &
       bind(c, name='amrex_eb_flux_reg_fineadd_dm')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, mlo, mhi, ratio, vlo, vhi, clo, chi
    integer, intent(in) :: nc
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in   ) :: dm(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3),nc)
    real(rt), intent(inout) :: d (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nc)
    real(rt)                :: cvol ( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    real(rt), intent(in   ) :: vfrac( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3))

    integer :: i,j,k,n, ii,jj,kk, ioff, joff, koff, iii, jjj, kkk
    real(rt) :: threshold

    threshold = reredistribution_threshold*(ratio(1)*ratio(2)*ratio(3))

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             cvol(i,j,k) = sum(vfrac(i*ratio(1):i*ratio(1)+ratio(1)-1,  &
                  &                  j*ratio(2):j*ratio(2)+ratio(2)-1, &
                  &                  k*ratio(3):k*ratio(3)+ratio(3)-1))
             if (cvol(i,j,k).gt.threshold) then
                cvol(i,j,k) = 1._rt/cvol(i,j,k)
             else
                cvol(i,j,k) = 0._rt
             end if
          end do
       end do
    end do

    do n = 1, nc
       do k = lo(3), hi(3)
          kk = k*ratio(3)
          do j = lo(2), hi(2)
             jj = j*ratio(2)

             do koff = 0, ratio(3)-1
                kkk = kk + koff
                do joff = 0, ratio(2)-1
                   jjj = jj + joff
                   do ioff = 0, ratio(1)-1

                      do i = lo(1), hi(1)
                         ii = i*ratio(1)
                         iii = ii + ioff
                         !$omp atomic
                         d(i,j,k,n) = d(i,j,k,n) + dm(iii,jjj,kkk,n)*cvol(i,j,k)
                         !$omp end atomic
                      end do

                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine amrex_eb_flux_reg_fineadd_dm

end module amrex_eb_flux_reg_3d_module
