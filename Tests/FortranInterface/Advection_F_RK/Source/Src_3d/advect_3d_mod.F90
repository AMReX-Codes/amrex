module advect_module

  use amrex_base_module

  implicit none
  private

  public :: advection_dudt

contains

subroutine advection_dudt(lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            dudt, du_lo, du_hi, &
     &            vx  , vx_lo, vx_hi, &
     &            vy  , vy_lo, vy_hi, &
     &            vz  , vz_lo, vz_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            flxz, fz_lo, fz_hi, &
     &            dx)

  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module, only : compute_flux_mol_3d

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  real(amrex_real), intent(in) :: dx(3)
  integer, intent(in) :: ui_lo(3), ui_hi(3)
  integer, intent(in) :: du_lo(3), du_hi(3)
  integer, intent(in) :: vx_lo(3), vx_hi(3)
  integer, intent(in) :: vy_lo(3), vy_hi(3)
  integer, intent(in) :: vz_lo(3), vz_hi(3)
  integer, intent(in) :: fx_lo(3), fx_hi(3)
  integer, intent(in) :: fy_lo(3), fy_hi(3)
  integer, intent(in) :: fz_lo(3), fz_hi(3)
  real(amrex_real), intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2),ui_lo(3):ui_hi(3))
  real(amrex_real), intent(  out) :: dudt(du_lo(1):du_hi(1),du_lo(2):du_hi(2),du_lo(3):du_hi(3))
  real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
  real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2),vy_lo(3):vy_hi(3))
  real(amrex_real), intent(in   ) :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
  real(amrex_real), intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
  real(amrex_real), intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
  real(amrex_real), intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))

  integer :: i, j, k
  integer :: glo(3), ghi(3)
  real(amrex_real) :: dxinv(3)

  ! Some compiler may not support 'contiguous'.  Remove it in that case.
  real(amrex_real), dimension(:,:,:), pointer, contiguous :: slope

  dxinv = 1.0_amrex_real/dx

  glo = lo - 1
  ghi = hi + 1

  call bl_allocate(slope,glo(1), ghi(1), glo(2), ghi(2), glo(3), ghi(3))

  ! We like to allocate these **pointers** here and then pass them to a function
  ! to remove their pointerness for performance, because normally pointers could
  ! be aliasing.  We need to use pointers instead of allocatable arrays because
  ! we like to use AMReX's bl_allocate to allocate memory instead of the intrinsic
  ! allocate.  Bl_allocate is much faster than allocate inside OMP.
  ! Note that one MUST CALL BL_DEALLOCATE.

  call compute_flux_mol_3d(lo, hi, &
                           uin, ui_lo, ui_hi, &
                           vx, vx_lo, vx_hi, &
                           vy, vy_lo, vy_hi, &
                           vz, vz_lo, vz_hi, &
                           flxx, fx_lo, fx_hi, &
                           flxy, fy_lo, fy_hi, &
                           flxz, fz_lo, fz_hi, &
                           slope, glo, ghi)

  ! Return dudt = -div(F) at the stage state.
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           dudt(i,j,k) = (flxx(i,j,k) - flxx(i+1,j,k)) * dxinv(1) &
                +         (flxy(i,j,k) - flxy(i,j+1,k)) * dxinv(2) &
                +         (flxz(i,j,k) - flxz(i,j,k+1)) * dxinv(3)
        enddo
     enddo
  enddo

  call bl_deallocate(slope)

end subroutine advection_dudt

end module advect_module
