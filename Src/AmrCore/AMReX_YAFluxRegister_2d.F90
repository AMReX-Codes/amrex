module amrex_ya_flux_reg_2d_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_ya_flux_reg_nd_module, only : crse_cell, crse_fine_boundary_cell, fine_cell
  implicit none
  private

  public :: amrex_ya_flux_reg_crseadd, amrex_ya_flux_reg_fineadd

contains

  subroutine amrex_ya_flux_reg_crseadd (lo, hi, d, dlo, dhi, flag, fglo, fghi, &
       fx, fxlo, fxhi, fy, fylo, fyhi, dx, dt, nc) &
       bind(c,name='amrex_ya_flux_reg_crseadd')
    integer, dimension(2), intent(in) :: lo, hi, dlo, dhi, fglo, fghi, fxlo, fxhi, fylo, fyhi
    integer, intent(in) :: nc
    integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2))
    real(rt), intent(in   ) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),nc)
    real(rt), intent(in   ) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),nc)
    real(rt), intent(inout) :: d ( dlo(1): dhi(1), dlo(2): dhi(2),nc)
    real(rt), intent(in) :: dx(2), dt

    integer :: i,j
    real(rt) :: dtdx, dtdy

    dtdx = dt/dx(1)
    dtdy = dt/dx(2)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (flag(i,j).eq.crse_fine_boundary_cell) then
             
             if (flag(i-1,j).eq.fine_cell) then
                d(i,j,:) = d(i,j,:) - dtdx*fx(i,j,:)
             else if (flag(i+1,j).eq.fine_cell) then
                d(i,j,:) = d(i,j,:) + dtdx*fx(i+1,j,:)
             end if
             
             if (flag(i,j-1).eq.fine_cell) then
                d(i,j,:) = d(i,j,:) - dtdy*fy(i,j,:)
             else if (flag(i,j+1).eq.fine_cell) then
                d(i,j,:) = d(i,j,:) + dtdy*fy(i,j+1,:)
             end if

          end if
       end do
    end do

  end subroutine amrex_ya_flux_reg_crseadd


  subroutine amrex_ya_flux_reg_fineadd (lo, hi, d, dlo, dhi, f, flo, fhi, dx, dt, nc, dir, side, ratio) &
       bind(c,name='amrex_ya_flux_reg_fineadd')
    integer, dimension(2), intent(in) :: lo, hi, dlo, dhi, flo, fhi, ratio
    integer, intent(in) :: nc, dir, side
    real(rt), intent(in) :: dx(2), dt
    real(rt), intent(in   ) :: f(flo(1):fhi(1),flo(2):fhi(2),nc)
    real(rt), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),nc)
    
    integer :: i,j,n,ii,jj,ioff,joff
    real(rt) :: fac

    ! dx is fine.
    ! lo and hi are also relative to the fine box.

    if (dir .eq. 0) then ! x-direction, lo(1) == hi(1)

       fac = dt / (dx(1)*(ratio(1)*ratio(2)))

       if (side .eq. 0) then ! lo-side

          i = lo(1)
          ii = (i+1)*ratio(1) !!!
          
          do n = 1, nc
             do j = lo(2), hi(2)
                do joff = 0, ratio(2)-1
                   jj = j*ratio(2)+joff
                   d(i,j,n) = d(i,j,n) - fac*f(ii,jj,n)
                end do
             end do
          end do

       else ! hi-side

          i = lo(1)
          ii = i*ratio(1) !!!

          do n = 1, nc
             do j = lo(2), hi(2)
                do joff = 0, ratio(2)-1
                   jj = j*ratio(2)+joff
                   d(i,j,n) = d(i,j,n) + fac*f(ii,jj,n)
                end do
             end do
          end do

       end if

    else  ! y-direction, lo(2) == hi(2)

       fac = dt / (dx(2)*(ratio(1)*ratio(2)))

       if (side .eq. 0) then ! lo-side

          j = lo(2)
          jj = (j+1)*ratio(2) !!!

          do n = 1, nc
             do ioff = 0, ratio(1)-1
                do i = lo(1), hi(1)
                   ii = i*ratio(1)+ioff
                   d(i,j,n) = d(i,j,n) - fac*f(ii,jj,n)
                end do
             end do
          end do
             
       else ! hi-side

          j = lo(2)
          jj = j*ratio(2) !!!

          do n = 1, nc
             do ioff = 0, ratio(1)-1
                do i = lo(1), hi(1)
                   ii = i*ratio(1)+ioff
                   d(i,j,n) = d(i,j,n) + fac*f(ii,jj,n)
                end do
             end do
          end do

       end if

    end if
  end subroutine amrex_ya_flux_reg_fineadd


end module amrex_ya_flux_reg_2d_module

