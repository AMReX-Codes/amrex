module compute_flux_module

  use amrex_base_module

  implicit none

  private

  public :: compute_flux_mol_3d

contains

  subroutine compute_flux_mol_3d(lo, hi, &
                                 phi, ph_lo, ph_hi, &
                                 umac, u_lo, u_hi, &
                                 vmac, v_lo, v_hi, &
                                 wmac, w_lo, w_hi, &
                                 flxx, fx_lo, fx_hi, &
                                 flxy, fy_lo, fy_hi, &
                                 flxz, fz_lo, fz_hi, &
                                 slope, glo, ghi)

    use slope_module, only: slopex, slopey, slopez

    integer, intent(in) :: lo(3), hi(3), glo(3), ghi(3)
    integer, intent(in) :: ph_lo(3), ph_hi(3)
    integer, intent(in) ::  u_lo(3),  u_hi(3)
    integer, intent(in) ::  v_lo(3),  v_hi(3)
    integer, intent(in) ::  w_lo(3),  w_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    real(amrex_real), intent(in   ) :: phi (ph_lo(1):ph_hi(1),ph_lo(2):ph_hi(2),ph_lo(3):ph_hi(3))
    real(amrex_real), intent(in   ) :: umac( u_lo(1): u_hi(1), u_lo(2): u_hi(2), u_lo(3): u_hi(3))
    real(amrex_real), intent(in   ) :: vmac( v_lo(1): v_hi(1), v_lo(2): v_hi(2), v_lo(3): v_hi(3))
    real(amrex_real), intent(in   ) :: wmac( w_lo(1): w_hi(1), w_lo(2): w_hi(2), w_lo(3): w_hi(3))
    real(amrex_real), intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3))
    real(amrex_real), intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3))
    real(amrex_real), intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3))
    real(amrex_real), dimension(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3)) :: slope

    integer :: i, j, k
    real(amrex_real) :: phiface

    call slopex(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             if (umac(i,j,k) .lt. 0.0_amrex_real) then
                phiface = phi(i,j,k) - 0.5_amrex_real*slope(i,j,k)
             else
                phiface = phi(i-1,j,k) + 0.5_amrex_real*slope(i-1,j,k)
             end if

             flxx(i,j,k) = umac(i,j,k)*phiface
          end do
       end do
    end do

    call slopey(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             if (vmac(i,j,k) .lt. 0.0_amrex_real) then
                phiface = phi(i,j,k) - 0.5_amrex_real*slope(i,j,k)
             else
                phiface = phi(i,j-1,k) + 0.5_amrex_real*slope(i,j-1,k)
             end if

             flxy(i,j,k) = vmac(i,j,k)*phiface
          end do
       end do
    end do

    call slopez(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    do       k = lo(3), hi(3)+1
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (wmac(i,j,k) .lt. 0.0_amrex_real) then
                phiface = phi(i,j,k) - 0.5_amrex_real*slope(i,j,k)
             else
                phiface = phi(i,j,k-1) + 0.5_amrex_real*slope(i,j,k-1)
             end if

             flxz(i,j,k) = wmac(i,j,k)*phiface
          end do
       end do
    end do

  end subroutine compute_flux_mol_3d

end module compute_flux_module
