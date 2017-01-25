
module warpx_laser_module

  use iso_c_binding
  use bl_fort_module, only : c_real
  use constants

  implicit none

contains

  subroutine warpx_gaussian_laser( np, Xp, Yp, t, &
      wavelength, waist, duration, t_peak, f, amplitude ) &
       bind(C, name="warpx_gaussian_laser")

    INTEGER(c_long), INTENT(IN) :: np
    REAL(c_real), INTENT(IN)    :: Xp(np),Yp(np)
    REAL(c_real), INTENT(IN)    :: t, t_peak, wavelength, duration, f, waist
    REAL(c_real), INTENT(INOUT) :: amplitude(np)

    INTEGER      :: i
    REAL(c_real) :: k0, oscillation_phase, temporal_exponent
    COMPLEX*16   :: diffract_factor, exp_argument, prefactor, &
                        inv_complex_waist_2, j=CMPLX(0., 1.)

    ! This function uses the complex expression of a Gaussian laser
    ! (Including Gouy phase and laser oscillations)

    ! Calculate a few factors which are independent of the macroparticle
    k0 = 2*pi/wavelength
    oscillation_phase = j * k0 * clight * ( t - t_peak )
    temporal_exponent = ( (t - t_peak) / duration )**2
    ! The coefficients below contain info about Gouy phase,
    ! laser diffraction, and phase front curvature
    diffract_factor = 1 - j * f * 2./(k0*waist**2)
    inv_complex_waist_2 = 1./( waist**2 * diffract_factor )

    ! Calculate amplitude prefactor
    prefactor = EXP( oscillation_phase - temporal_exponent )
    ! Because diffract_factor is a complex, the code below takes into
    ! account the impact of the dimensionality on both the Gouy phase
    ! and the amplitude of the laser
#if (BL_SPACEDIM == 3)
    prefactor = prefactor / diffract_factor
#elif (BL_SPACEDIM == 2)
    prefactor = prefactor / POW( diffract_factor, 0.5 )
#endif

    ! Loop through the macroparticle to calculate the proper amplitude
    DO i = 1, np
      exp_argument = - ( Xp(i)*Xp(i) + Yp(i)*Yp(i) ) * inv_complex_waist_2
      amplitude(i) = REALPART( prefactor * EXP( exp_argument ) )
    ENDDO

  end subroutine warpx_gaussian_laser

end module warpx_laser_module
