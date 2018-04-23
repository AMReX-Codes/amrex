
module warpx_laser_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use constants
  use parser_wrapper

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

  ! Harris function for the laser temporal profile
  subroutine warpx_harris_laser( np, Xp, Yp, t, &
      wavelength, e_max, waist, duration, f, amplitude ) &
       bind(C, name="warpx_harris_laser")

       integer(c_long), intent(in) :: np
       real(amrex_real), intent(in)    :: Xp(np),Yp(np)
       real(amrex_real), intent(in)    :: e_max, t, wavelength, duration, f, waist
       real(amrex_real), intent(inout) :: amplitude(np)

       integer(c_long)  :: i
       real(amrex_real) :: space_envelope, time_envelope, arg_osc, arg_env
       real(amrex_real) :: omega0, zR, wz, inv_Rz, oscillations, inv_wz_2

    ! This function uses the Harris function as the temporal profile of the pulse
    omega0 = 2*pi*clight/wavelength
    zR = pi * waist**2 / wavelength
    wz = waist * sqrt(1. + f**2/zR**2)
    inv_wz_2 = 1./wz**2
    if (f == 0.) then
      inv_Rz = 0.
    else
      inv_Rz = -f / ( f**2 + zR**2 )
    end if

    arg_env = 2.*pi*t/duration

    ! time envelope is given by the Harris function
    time_envelope = 0.
    if (t < duration) then
      time_envelope = 1./32. * (10. - 15.*cos(arg_env) + 6.*cos(2.*arg_env) - cos(3.*arg_env))
    else
      time_envelope = 0.
    end if

    ! Loop through the macroparticle to calculate the proper amplitude
    do i = 1, np
      space_envelope = exp(- ( Xp(i)*Xp(i) + Yp(i)*Yp(i) ) * inv_wz_2)
      arg_osc = omega0*t - omega0/clight*(Xp(i)*Xp(i) + Yp(i)*Yp(i)) * inv_Rz / 2
      oscillations = cos(arg_osc)
      amplitude(i) = e_max * time_envelope * space_envelope * oscillations
    enddo

  end subroutine warpx_harris_laser

  ! Parse function from the input script for the laser temporal profile
  subroutine parse_function_laser( np, Xp, Yp, t, e_max, amplitude, parser_instance_number ) bind(C, name="parse_function_laser")
       integer(c_long), intent(in) :: np
       real(amrex_real), intent(in)    :: Xp(np),Yp(np)
       real(amrex_real), intent(in)    :: e_max, t
       real(amrex_real), intent(inout) :: amplitude(np)
       INTEGER, value, INTENT(IN) :: parser_instance_number
       integer(c_long)  :: i
       INTEGER, PARAMETER :: nvar_parser = 3
       REAL(amrex_real) :: list_var(1:nvar_parser)
    ! Loop through the macroparticle to calculate the proper amplitude
    do i = 1, np
      list_var = [Xp(i), Yp(i), t]
      amplitude(i) = e_max * parser_evaluate_function(list_var, nvar_parser, parser_instance_number)
    enddo
  end subroutine parse_function_laser

end module warpx_laser_module
