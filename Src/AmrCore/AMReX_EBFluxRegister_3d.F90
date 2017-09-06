module amrex_eb_flux_reg_3d_module

  implicit none
  private

  integer, parameter :: crse_cell = 0
  integer, parameter :: crse_fine_boundary_cell = 1
  integer, parameter :: fine_cell = 2

  public :: amrex_eb_flux_reg_crseadd

contains

  subroutine amrex_eb_flux_reg_crseadd (lo, hi, d, dlo, dhi, flag, fglo, fghi, &
       fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, dx, dt, nc) &
       bind(c,name='amrex_eb_flux_reg_crseadd')
    use amrex_fort_module, only : rt => amrex_real
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
                   d(i,j,k,:) = d(i,j,k,:) + dtdz*fz(i,j,k,:)
                end if

             end if
          end do
       end do
    end do

  end subroutine amrex_eb_flux_reg_crseadd

end module amrex_eb_flux_reg_3d_module
