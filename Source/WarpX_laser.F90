
module warpx_laser_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use constants

  implicit none

contains

  subroutine warpx_gaussian_laser( np, Xp, Yp, t, &
      wavelength, e_max, waist, duration, t_peak, f, amplitude ) &
       bind(C, name="warpx_gaussian_laser")

    integer(c_long), intent(in) :: np
    real(amrex_real), intent(in)    :: Xp(np),Yp(np)
    real(amrex_real), intent(in)    :: e_max, t, t_peak, wavelength, duration, f, waist
    real(amrex_real), intent(inout) :: amplitude(np)

    integer(c_long)  :: i
    real(amrex_real) :: k0, oscillation_phase, temporal_exponent
    complex*16       :: diffract_factor, exp_argument, prefactor, &
                        inv_complex_waist_2
    complex*16, parameter :: j=cmplx(0., 1.)

    ! This function uses the complex expression of a Gaussian laser
    ! (Including Gouy phase and laser oscillations)

    ! Calculate a few factors which are independent of the macroparticle
    k0 = 2*pi/wavelength
    oscillation_phase = k0 * clight * ( t - t_peak )
    temporal_exponent = ( (t - t_peak) / duration )**2
    ! The coefficients below contain info about Gouy phase,
    ! laser diffraction, and phase front curvature
    diffract_factor = 1 + j * f * 2./(k0*waist**2)
    inv_complex_waist_2 = 1./( waist**2 * diffract_factor )

    ! Calculate amplitude prefactor
    prefactor = e_max * exp( j * oscillation_phase - temporal_exponent )

    ! Because diffract_factor is a complex, the code below takes into
    ! account the impact of the dimensionality on both the Gouy phase
    ! and the amplitude of the laser
#if (BL_SPACEDIM == 3)
    prefactor = prefactor / diffract_factor
#elif (BL_SPACEDIM == 2)
    prefactor = prefactor / sqrt(diffract_factor)
#endif

    ! Loop through the macroparticle to calculate the proper amplitude
    do i = 1, np
      exp_argument = - ( Xp(i)*Xp(i) + Yp(i)*Yp(i) ) * inv_complex_waist_2
      amplitude(i) = DREAL( prefactor * exp( exp_argument ) )
    enddo

  end subroutine warpx_gaussian_laser

end module warpx_laser_module
