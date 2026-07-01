module compute_flux_module

  use amrex_base_module

  implicit none

  private

  public :: compute_flux_mol_2d

contains

  subroutine compute_flux_mol_2d(lo, hi, &
                                 phi, ph_lo, ph_hi, &
                                 umac, u_lo, u_hi, &
                                 vmac, v_lo, v_hi, &
                                 flxx, fx_lo, fx_hi, &
                                 flxy, fy_lo, fy_hi, &
                                 slope, glo, ghi)

    use slope_module, only: slopex, slopey

    integer, intent(in) :: lo(2), hi(2), glo(2), ghi(2)
    integer, intent(in) :: ph_lo(2), ph_hi(2)
    integer, intent(in) ::  u_lo(2),  u_hi(2)
    integer, intent(in) ::  v_lo(2),  v_hi(2)
    integer, intent(in) :: fx_lo(2), fx_hi(2)
    integer, intent(in) :: fy_lo(2), fy_hi(2)
    real(amrex_real), intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2))
    real(amrex_real), intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2))
    real(amrex_real), intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2))
    real(amrex_real), intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
    real(amrex_real), intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))
    real(amrex_real), dimension(glo(1):ghi(1),glo(2):ghi(2)) :: slope

    integer :: i, j
    real(amrex_real) :: phiface

    call slopex(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          if (umac(i,j) .lt. 0.0_amrex_real) then
             phiface = phi(i,j) - 0.5_amrex_real*slope(i,j)
          else
             phiface = phi(i-1,j) + 0.5_amrex_real*slope(i-1,j)
          end if

          flxx(i,j) = phiface*umac(i,j)
       end do
    end do

    call slopey(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          if (vmac(i,j) .lt. 0.0_amrex_real) then
             phiface = phi(i,j) - 0.5_amrex_real*slope(i,j)
          else
             phiface = phi(i,j-1) + 0.5_amrex_real*slope(i,j-1)
          end if

          flxy(i,j) = phiface*vmac(i,j)
       end do
    end do

  end subroutine compute_flux_mol_2d

end module compute_flux_module
