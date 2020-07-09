module compute_flux_module

  use amrex_base_module

  implicit none

  private

  public :: compute_flux

contains

  subroutine compute_flux(lo, hi, &
                          phi,ph_lo,ph_hi, &
                          umac,  u_lo,  u_hi, &
                          vmac,  v_lo,  v_hi, &
                          flxx, fx_lo, fx_hi, &
                          flxy, fy_lo, fy_hi, &
                          dt, dx)

    use slope_module, only: slopex, slopey

    integer, intent(in) :: lo(2), hi(2)
    real(amrex_real), intent(in) :: dt, dx(2)
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
         
    integer :: i, j, k, glo(2), ghi(2)
    real(amrex_real) :: hdtdx(2), umax, vmax
    real(amrex_real), dimension(:,:), pointer, contiguous :: phix_1d, phiy_1d, phix, phiy, slope    

    ! check if CFL condition is violated.
    umax = maxval(abs(umac))
    vmax = maxval(abs(vmac))
    if ( umax*dt .ge. dx(1) .or. &
         vmax*dt .ge. dx(2) ) then
       print *, "umax = ", umax, ", vmax = ", vmax, ", dt = ", dt, ", dx = ", dx
       call amrex_error("CFL violation. Use smaller adv.cfl.")
    end if

    hdtdx = 0.5_amrex_real*(dt/dx)

    glo = lo - 1
    ghi = hi + 1

    ! edge states
    call amrex_allocate(phix_1d, glo, ghi)
    call amrex_allocate(phiy_1d, glo, ghi)
    call amrex_allocate(phix   , glo, ghi)
    call amrex_allocate(phiy   , glo, ghi)
    ! slope
    call amrex_allocate(slope  , glo, ghi)

    call slopex(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute phi on x faces using umac to upwind; ignore transverse terms
    do    j = lo(2)-1, hi(2)+1
       do i = lo(1)  , hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             phix_1d(i,j) = phi(i  ,j) - (0.5d0 + hdtdx(1)*umac(i,j))*slope(i  ,j)
          else
             phix_1d(i,j) = phi(i-1,j) + (0.5d0 - hdtdx(1)*umac(i,j))*slope(i-1,j)
          end if

       end do
    end do

    call slopey(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)

    ! compute phi on y faces using umac to upwind; ignore transverse terms
    do    j = lo(2)  , hi(2)+1
       do i = lo(1)-1, hi(1)+1

          if (vmac(i,j) .lt. 0.d0) then
             phiy_1d(i,j) = phi(i,j  ) - (0.5d0 + hdtdx(2)*vmac(i,j))*slope(i,j  )
          else
             phiy_1d(i,j) = phi(i,j-1) + (0.5d0 - hdtdx(2)*vmac(i,j))*slope(i,j-1)
          end if

       end do
    end do

    ! update phi on x faces by adding in y-transverse terms
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1

          if (umac(i,j) .lt. 0.d0) then
             phix(i,j) = phix_1d(i,j) &
                  - hdtdx(2)*( 0.5d0*(vmac(i  ,j+1)+vmac(i  ,j)) * (phiy_1d(i  ,j+1)-phiy_1d(i  ,j)) )
          else
             phix(i,j) = phix_1d(i,j) &
                  - hdtdx(2)*( 0.5d0*(vmac(i-1,j+1)+vmac(i-1,j)) * (phiy_1d(i-1,j+1)-phiy_1d(i-1,j)) )
          end if

          ! compute final x-fluxes
          flxx(i,j) = phix(i,j)*umac(i,j)

       end do
    end do

    ! update phi on y faces by adding in x-transverse terms
    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)

          if (vmac(i,j) .lt. 0.d0) then
             phiy(i,j) = phiy_1d(i,j) &
                  - hdtdx(1)*( 0.5d0*(umac(i+1,j  )+umac(i,j  )) * (phix_1d(i+1,j  )-phix_1d(i,j  )) )
          else
             phiy(i,j) = phiy_1d(i,j) &
                  - hdtdx(1)*( 0.5d0*(umac(i+1,j-1)+umac(i,j-1)) * (phix_1d(i+1,j-1)-phix_1d(i,j-1)) )
          end if

          ! compute final y-fluxes
          flxy(i,j) = phiy(i,j)*vmac(i,j)

       end do
    end do

    ! Final fluxes
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)+1
          flxx(i,j) = phix(i,j) * umac(i,j)
       end do
    end do
    !
    do    j = lo(2), hi(2)+1
       do i = lo(1), hi(1)
          flxy(i,j) = phiy(i,j) * vmac(i,j)
       end do
    end do

    call amrex_deallocate(phix_1d)
    call amrex_deallocate(phiy_1d)
    call amrex_deallocate(phix)
    call amrex_deallocate(phiy)
    call amrex_deallocate(slope)
    
  end subroutine compute_flux

end module compute_flux_module
