module amrex_ya_flux_reg_1d_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_ya_flux_reg_nd_module, only : crse_cell, crse_fine_boundary_cell, fine_cell
  implicit none
  private

  public :: amrex_ya_flux_reg_crseadd, amrex_ya_flux_reg_fineadd

contains

  subroutine amrex_ya_flux_reg_crseadd (lo, hi, d, dlo, dhi, flag, fglo, fghi, &
       fx, fxlo, fxhi, dx, dt, nc) &
       bind(c,name='amrex_ya_flux_reg_crseadd')
    integer, dimension(1), intent(in) :: lo, hi, dlo, dhi, fglo, fghi, fxlo, fxhi
    integer, intent(in) :: nc
    integer, intent(in) :: flag(fglo(1):fghi(1))
    real(rt), intent(in   ) :: fx(fxlo(1):fxhi(1),nc)
    real(rt), intent(inout) :: d ( dlo(1): dhi(1),nc)
    real(rt), intent(in) :: dx(1), dt

    integer :: i
    real(rt) :: dtdx

    dtdx = dt/dx(1)

    do i = lo(1), hi(1)
       if (flag(i).eq.crse_fine_boundary_cell) then
          
          if (flag(i-1).eq.fine_cell) then
             d(i,:) = d(i,:) - dtdx*fx(i,:)
          else if (flag(i+1).eq.fine_cell) then
             d(i,:) = d(i,:) + dtdx*fx(i+1,:)
          end if

       end if
    end do

  end subroutine amrex_ya_flux_reg_crseadd


  subroutine amrex_ya_flux_reg_fineadd (lo, hi, d, dlo, dhi, f, flo, fhi, dx, dt, nc, dir, side, ratio) &
       bind(c,name='amrex_ya_flux_reg_fineadd')
    integer, dimension(1), intent(in) :: lo, hi, dlo, dhi, flo, fhi, ratio
    integer, intent(in) :: nc, dir, side
    real(rt), intent(in) :: dx(1), dt
    real(rt), intent(in   ) :: f(flo(1):fhi(1),nc)
    real(rt), intent(inout) :: d(dlo(1):dhi(1),nc)
    
    integer :: i,n,ii,ioff
    real(rt) :: fac

    ! dx is fine.
    ! lo and hi are also relative to the fine box.

    ! lo(1) == hi(1)

    fac = dt / (dx(1)*ratio(1))

    if (side .eq. 0) then ! lo-side

       i = lo(1)
       ii = (i+1)*ratio(1) !!!          
       d(i,:) = d(i,:) - fac*f(ii,:)

    else ! hi-side

       i = lo(1)
       ii = i*ratio(1) !!!
       d(i,:) = d(i,:) + fac*f(ii,:)

    end if

  end subroutine amrex_ya_flux_reg_fineadd


end module amrex_ya_flux_reg_1d_module

