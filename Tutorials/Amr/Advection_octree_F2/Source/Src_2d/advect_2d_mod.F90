
module advect_module

  use amrex_base_module

  implicit none
  private
  
  public :: advect

contains

  subroutine advect(lo, hi, &
       &            uin , ui_lo, ui_hi, &
       &            uout, uo_lo, uo_hi, &
       &            flxx, fx_lo, fx_hi, &
       &            flxy, fy_lo, fy_hi, &
       &            dx,dt)

    use compute_flux_module, only : compute_flux

    integer, intent(in) :: lo(2), hi(2)
    real(amrex_real), intent(in) :: dx(2), dt
    integer, intent(in) :: ui_lo(2), ui_hi(2)
    integer, intent(in) :: uo_lo(2), uo_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    real(amrex_real), intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2))
    real(amrex_real), intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
    real(amrex_real), intent(in   ) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    real(amrex_real), intent(in   ) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    
    integer :: i, j
    integer :: glo(2), ghi(2)
    real(amrex_real) :: dtdx(2), umax, vmax

    dtdx = dt/dx

    ! Do a conservative update
    do    j = lo(2),hi(2)
       do i = lo(1),hi(1)
          uout(i,j) = uin(i,j) + &
               ( (flxx(i,j) - flxx(i+1,j)) * dtdx(1) &
               + (flxy(i,j) - flxy(i,j+1)) * dtdx(2) )
       enddo
    enddo
    
  end subroutine advect

end module advect_module
