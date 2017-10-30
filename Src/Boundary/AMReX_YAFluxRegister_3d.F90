module amrex_ya_flux_reg_3d_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_ya_flux_reg_nd_module, only : crse_cell, crse_fine_boundary_cell, fine_cell
  implicit none
  private

  public :: amrex_ya_flux_reg_crseadd, amrex_ya_flux_reg_fineadd

contains

  subroutine amrex_ya_flux_reg_crseadd (lo, hi, d, dlo, dhi, flag, fglo, fghi, &
       fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, dx, dt, nc) &
       bind(c,name='amrex_ya_flux_reg_crseadd')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, fglo, fghi, fxlo, fxhi, fylo, fyhi, fzlo, fzhi
    integer, intent(in) :: nc
    integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    real(rt), intent(in   ) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),nc)
    real(rt), intent(in   ) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),nc)
    real(rt), intent(in   ) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),nc)
    real(rt), intent(inout) :: d ( dlo(1): dhi(1), dlo(2): dhi(2), dlo(3): dhi(3),nc)
    real(rt), intent(in) :: dx(3), dt

    integer :: i,j,k
    real(rt) :: dtdx, dtdy, dtdz

    dtdx = dt/dx(1)
    dtdy = dt/dx(2)
    dtdz = dt/dx(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (flag(i,j,k).eq.crse_fine_boundary_cell) then

                if (flag(i-1,j,k).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) - dtdx*fx(i,j,k,:)
                else if (flag(i+1,j,k).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) + dtdx*fx(i+1,j,k,:)
                end if

                if (flag(i,j-1,k).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) - dtdy*fy(i,j,k,:)
                else if (flag(i,j+1,k).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) + dtdy*fy(i,j+1,k,:)
                end if

                if (flag(i,j,k-1).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) - dtdz*fz(i,j,k,:)
                else if (flag(i,j,k+1).eq.fine_cell) then
                   d(i,j,k,:) = d(i,j,k,:) + dtdz*fz(i,j,k+1,:)
                end if

             end if
          end do
       end do
    end do

  end subroutine amrex_ya_flux_reg_crseadd


  subroutine amrex_ya_flux_reg_fineadd (lo, hi, d, dlo, dhi, f, flo, fhi, dx, dt, nc, dir, side, ratio) &
       bind(c,name='amrex_ya_flux_reg_fineadd')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, flo, fhi, ratio
    integer, intent(in) :: nc, dir, side
    real(rt), intent(in) :: dx(3), dt
    real(rt), intent(in   ) :: f(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),nc)
    real(rt), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nc)
    
    integer :: i,j,k,n,ii,jj,kk,ioff,joff,koff
    real(rt) :: fac

    ! dx is fine.
    ! lo and hi are also relative to the fine box.

    if (dir .eq. 0) then ! x-direction, lo(1) == hi(1)

       fac = dt / (dx(1)*(ratio(1)*ratio(2)*ratio(3)))

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
                         !$omp atomic
                         d(i,j,k,n) = d(i,j,k,n) - fac*f(ii,jj,kk,n)
                         !$omp end atomic
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
                         !$omp atomic
                         d(i,j,k,n) = d(i,j,k,n) + fac*f(ii,jj,kk,n)
                         !$omp end atomic
                      end do
                   end do
                end do
             end do
          end do

       end if

    else if (dir .eq. 1) then ! y-direction, lo(2) == hi(2)

       fac = dt / (dx(2)*(ratio(1)*ratio(2)*ratio(3)))

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
                         !$omp atomic
                         d(i,j,k,n) = d(i,j,k,n) - fac*f(ii,jj,kk,n)
                         !$omp end atomic
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
                         !$omp atomic
                         d(i,j,k,n) = d(i,j,k,n) + fac*f(ii,jj,kk,n)
                         !$omp end atomic
                      end do
                   end do
                end do
             end do
          end do

       end if

    else  ! z-direction, lo(3) == hi(3)

       fac = dt / (dx(3)*(ratio(1)*ratio(2)*ratio(3)))

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
                         !$omp atomic
                         d(i,j,k,n) = d(i,j,k,n) - fac*f(ii,jj,kk,n)
                         !$omp end atomic
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
                         !$omp atomic
                         d(i,j,k,n) = d(i,j,k,n) + fac*f(ii,jj,kk,n)
                         !$omp end atomic
                      end do
                   end do
                end do
             end do
          end do

       end if
    end if
  end subroutine amrex_ya_flux_reg_fineadd


end module amrex_ya_flux_reg_3d_module
