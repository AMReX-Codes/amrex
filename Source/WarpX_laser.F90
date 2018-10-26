#if AMREX_USE_GPU && __PGIC__ <= 18 && __PGIC_MINOR__ <= 7
#define USE_ACC_LASER_SUBROUTINE
! Define a separate subroutine that contains the openacc pragmas
! This is because fortran subroutines containing openacc pramas
! ignore the `bind(C,` statement, and so cannot directly be called from C++
! (at least until pgi version 18.7)
#endif

module warpx_laser_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use constants, only : clight, pi
  use parser_wrapper, only : parser_evaluate_function

  implicit none

contains

  subroutine warpx_gaussian_laser( np, Xp, Yp, t, &
      wavelength, e_max, waist, duration, t_peak, f, amplitude, &
      zeta, beta, phi2 ) bind(C, name="warpx_gaussian_laser")

    integer(c_long), intent(in) :: np
    real(amrex_real), intent(in)    :: Xp(np),Yp(np)
    real(amrex_real), intent(in)    :: e_max, t, t_peak, wavelength, duration, f, waist
    real(amrex_real), intent(in)    :: zeta, beta, phi2
    real(amrex_real), intent(inout) :: amplitude(np)

    integer(c_long)  :: i
    real(amrex_real) :: k0, oscillation_phase, inv_tau2
    complex*16       :: diffract_factor, exp_argument, prefactor, &
                        inv_complex_waist_2, stretch_factor, &
                        stc_exponent, stcfactor
    complex*16, parameter :: j=cmplx(0., 1.)

#ifdef USE_ACC_LASER_SUBROUTINE
    call warpx_gaussian_laser_acc( np, Xp, Yp, t, &
      wavelength, e_max, waist, duration, t_peak, f, amplitude, &
      zeta, beta, phi2 )

  end subroutine warpx_gaussian_laser
    
  subroutine warpx_gaussian_laser_acc( np, Xp, Yp, t, &
      wavelength, e_max, waist, duration, t_peak, f, amplitude, &
      zeta, beta, phi2 )

    integer(c_long), intent(in) :: np
    real(amrex_real), intent(in)    :: Xp(np),Yp(np)
    real(amrex_real), intent(in)    :: e_max, t, t_peak, wavelength, duration, f, waist
    real(amrex_real), intent(in)    :: zeta, beta, phi2
    real(amrex_real), intent(inout) :: amplitude(np)

    integer(c_long)  :: i
    real(amrex_real) :: k0, oscillation_phase, inv_tau2
    complex*16       :: diffract_factor, exp_argument, prefactor, &
                        inv_complex_waist_2, stretch_factor, &
                        stc_exponent, stcfactor
    complex*16, parameter :: j=cmplx(0., 1.)
#endif
    
    ! This function uses the complex expression of a Gaussian laser
    ! (Including Gouy phase and laser oscillations)

    ! Calculate a few factors which are independent of the macroparticle
    k0 = 2*pi/wavelength
    inv_tau2 = 1. / duration**2
    oscillation_phase = k0 * clight * ( t - t_peak )
    ! The coefficients below contain info about Gouy phase,
    ! laser diffraction, and phase front curvature
    diffract_factor = 1 + j * f * 2./(k0*waist**2)
    inv_complex_waist_2 = 1./( waist**2 * diffract_factor )
    
    ! Time stretching due to STCs and phi2 complex envelope
    ! (1 if zeta=0, beta=0, phi2=0) 
   stretch_factor = 1. + &
        4*(zeta + beta*f)**2 * (inv_tau2*inv_complex_waist_2) + &
        2*j*(phi2 - beta**2*k0*f) * inv_tau2
    
    ! Amplitude and monochromatic oscillations
    prefactor = e_max * exp( j * oscillation_phase )

    ! Because diffract_factor is a complex, the code below takes into
    ! account the impact of the dimensionality on both the Gouy phase
    ! and the amplitude of the laser
#if (AMREX_SPACEDIM == 3)
    prefactor = prefactor / diffract_factor
#elif (AMREX_SPACEDIM == 2)
    prefactor = prefactor / sqrt(diffract_factor)
#endif

    ! Loop through the macroparticle to calculate the proper amplitude
    !$acc parallel deviceptr(amplitude, Xp, Yp)
    !$acc loop
    do i = 1, np
      ! Exp argument for the temporal gaussian envelope + STCs
      stc_exponent = 1./stretch_factor * inv_tau2 * &
          (t - t_peak - beta*k0*Xp(i) - 2*j*Xp(i)*( zeta - beta*f ) * &
          inv_complex_waist_2)**2 
      ! stcfactor = everything but complex transverse envelope 
      stcfactor = prefactor * exp( - stc_exponent )
      ! Exp argument for transverse envelope
      exp_argument = - ( Xp(i)*Xp(i) + Yp(i)*Yp(i) ) * inv_complex_waist_2
      ! stcfactor + transverse envelope
      amplitude(i) = DREAL( stcfactor * exp( exp_argument ) )
   enddo
   !$acc end loop
   !$acc end parallel

#ifdef USE_ACC_LASER_SUBROUTINE
  end subroutine warpx_gaussian_laser_acc
#else
  end subroutine warpx_gaussian_laser
#endif
  
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

#ifdef USE_ACC_LASER_SUBROUTINE
       call warpx_harris_laser_acc( np, Xp, Yp, t, &
            wavelength, e_max, waist, duration, f, amplitude )

  end subroutine warpx_harris_laser
       
  subroutine warpx_harris_laser_acc( np, Xp, Yp, t, &
          wavelength, e_max, waist, duration, f, amplitude )

       integer(c_long), intent(in) :: np
       real(amrex_real), intent(in)    :: Xp(np),Yp(np)
       real(amrex_real), intent(in)    :: e_max, t, wavelength, duration, f, waist
       real(amrex_real), intent(inout) :: amplitude(np)

       integer(c_long)  :: i
       real(amrex_real) :: space_envelope, time_envelope, arg_osc, arg_env
       real(amrex_real) :: omega0, zR, wz, inv_Rz, oscillations, inv_wz_2
#endif

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
    !$acc parallel deviceptr(amplitude, Xp, Yp)
    !$acc loop
    do i = 1, np
      space_envelope = exp(- ( Xp(i)*Xp(i) + Yp(i)*Yp(i) ) * inv_wz_2)
      arg_osc = omega0*t - omega0/clight*(Xp(i)*Xp(i) + Yp(i)*Yp(i)) * inv_Rz / 2
      oscillations = cos(arg_osc)
      amplitude(i) = e_max * time_envelope * space_envelope * oscillations
    enddo
    !$acc end loop
    !$acc end parallel

#ifdef USE_ACC_LASER_SUBROUTINE
  end subroutine warpx_harris_laser_acc
#else
  end subroutine warpx_harris_laser
#endif
  
  ! Parse function from the input script for the laser temporal profile
  subroutine parse_function_laser( np, Xp, Yp, t, amplitude, parser_instance_number ) bind(C, name="parse_function_laser")
       integer(c_long), intent(in) :: np
       real(amrex_real), intent(in)    :: Xp(np),Yp(np)
       real(amrex_real), intent(in)    :: t
       real(amrex_real), intent(inout) :: amplitude(np)
       INTEGER, value, INTENT(IN) :: parser_instance_number
       integer(c_long)  :: i
       INTEGER, PARAMETER :: nvar_parser = 3
       REAL(amrex_real) :: list_var(1:nvar_parser)
    ! Loop through the macroparticle to calculate the proper amplitude
    do i = 1, np
      list_var = [Xp(i), Yp(i), t]
      amplitude(i) = parser_evaluate_function(list_var, nvar_parser, parser_instance_number)
    enddo
  end subroutine parse_function_laser


  subroutine calculate_laser_plane_coordinates(np, xp, yp, zp, &
       plane_Xp, plane_Yp, u_Xx, u_Xy, u_Xz, u_Yx, u_Yy, u_Yz, &
       positionx, positiony, positionz ) &
       bind(C, name="calculate_laser_plane_coordinates")
    integer(c_long), intent(in) :: np
    real(amrex_real), intent(in) :: xp(np), yp(np), zp(np)
    real(amrex_real), intent(in) :: u_Xx, u_Xy, u_Xz
    real(amrex_real), intent(in) :: u_Yx, u_Yy, u_Yz
    real(amrex_real), intent(in) :: positionx, positiony, positionz    
    real(amrex_real), intent(out) :: plane_Xp(np), plane_Yp(np)
    integer :: ip

#ifdef USE_ACC_LASER_SUBROUTINE
    call calculate_laser_plane_coordinates_acc(np, xp, yp, zp, &
       plane_Xp, plane_Yp, u_Xx, u_Xy, u_Xz, u_Yx, u_Yy, u_Yz, &
       positionx, positiony, positionz )

  end subroutine calculate_laser_plane_coordinates

   subroutine calculate_laser_plane_coordinates_acc(np, xp, yp, zp, &
       plane_Xp, plane_Yp, u_Xx, u_Xy, u_Xz, u_Yx, u_Yy, u_Yz, &
       positionx, positiony, positionz )
    integer(c_long), intent(in) :: np
    real(amrex_real), intent(in) :: xp(np), yp(np), zp(np)
    real(amrex_real), intent(in) :: u_Xx, u_Xy, u_Xz
    real(amrex_real), intent(in) :: u_Yx, u_Yy, u_Yz
    real(amrex_real), intent(in) :: positionx, positiony, positionz    
    real(amrex_real), intent(out) :: plane_Xp(np), plane_Yp(np)
    integer :: ip
#endif

    !$acc parallel deviceptr(xp, yp, zp, plane_Xp, plane_Yp)
    !$acc loop
    do ip = 1, np
#if (BL_SPACEDIM == 3)
       plane_Xp(ip) = u_Xx*(xp(ip) - positionx) &
                   + u_Xy*(yp(ip) - positiony) &
                   + u_Xz*(zp(ip) - positionz)
       plane_Yp(ip) = u_Yx*(xp(ip) - positionx) &
                   + u_Yy*(yp(ip) - positiony) &
                   + u_Yz*(zp(ip) - positionz)
   #elif (AMREX_SPACEDIM == 2)
       plane_Xp(ip) = u_Xx*(xp(ip) - positionx) &
                   + u_Xz*(zp(ip) - positionz)
       plane_Yp(ip) = 0
#endif
    enddo
    !$acc end loop
    !$acc end parallel

#ifdef USE_ACC_LASER_SUBROUTINE
  end subroutine calculate_laser_plane_coordinates_acc
#else
  end subroutine calculate_laser_plane_coordinates
#endif

  subroutine update_laser_particle( np, &
    xp, yp, zp, uxp, uyp, uzp, giv, wp, amplitude_E, p_Xx, p_Xy, p_Xz, &
    nvecx, nvecy, nvecz, mobility, dt, c, beta_boost, gamma_boost ) &
    bind(C, name="update_laser_particle")

    integer(c_long), intent(in) :: np
    real(amrex_real), intent(in) :: amplitude_E(np), wp(np)
    real(amrex_real), intent(inout) :: xp(np), yp(np), zp(np)
    real(amrex_real), intent(out) :: uxp(np), uyp(np), uzp(np), giv(np)
    real(amrex_real), intent(in) :: p_Xx, p_Xy, p_Xz
    real(amrex_real), intent(in) :: nvecx, nvecy, nvecz
    real(amrex_real), intent(in) :: c, beta_boost, gamma_boost, mobility, dt

#ifdef USE_ACC_LASER_SUBROUTINE
    call update_laser_particle_acc( np, &
    xp, yp, zp, uxp, uyp, uzp, giv, wp, amplitude_E, p_Xx, p_Xy, p_Xz, &
    nvecx, nvecy, nvecz, mobility, dt, c, beta_boost, gamma_boost )

  end subroutine update_laser_particle

  subroutine update_laser_particle_acc( np, &
    xp, yp, zp, uxp, uyp, uzp, giv, wp, amplitude_E, p_Xx, p_Xy, p_Xz, &
    nvecx, nvecy, nvecz, mobility, dt, c, beta_boost, gamma_boost )

    integer(c_long), intent(in) :: np
    real(amrex_real), intent(in) :: amplitude_E(np), wp(np)
    real(amrex_real), intent(inout) :: xp(np), yp(np), zp(np)
    real(amrex_real), intent(out) :: uxp(np), uyp(np), uzp(np), giv(np)    
    real(amrex_real), intent(in) :: p_Xx, p_Xy, p_Xz
    real(amrex_real), intent(in) :: nvecx, nvecy, nvecz
    real(amrex_real), intent(in) :: c, beta_boost, gamma_boost, mobility, dt
#endif
    
    real(amrex_real)     :: vx, vy, vz, v_over_c, sign_charge, gamma
    integer  :: ip

    !$acc parallel deviceptr(amplitude_E, wp, xp, yp, zp, uxp, uyp, uzp, giv)
    !$acc loop
    do ip = 1, np
       ! Calculate the velocity according to the amplitude of E       
       sign_charge = SIGN( 1.0, wp(ip) )
       v_over_c = sign_charge * mobility * amplitude_E(ip)
       ! The velocity is along the laser polarization p_X
       vx = c * v_over_c * p_Xx
       vy = c * v_over_c * p_Xy
       vz = c * v_over_c * p_Xz
       ! When running in the boosted-frame, their is additional velocity along nvec
       if (gamma_boost > 1.) then
          vx = vx - c * beta_boost * nvecx
          vy = vy - c * beta_boost * nvecy
          vz = vz - c * beta_boost * nvecz
       endif
       ! Get the corresponding momenta
       gamma = gamma_boost/sqrt(1 - v_over_c*v_over_c)
       giv(ip) = 1./gamma
       uxp(ip) = gamma * vx
       uyp(ip) = gamma * vy
       uzp(ip) = gamma * vz
       ! Push the the particle positions
       xp(ip) = xp(ip) + vx * dt
#if (BL_SPACEDIM == 3)
       yp(ip) = yp(ip) + vy * dt
#endif
       zp(ip) = zp(ip) + vz * dt
    enddo
    !$acc end loop
    !$acc end parallel

#ifdef USE_ACC_LASER_SUBROUTINE
  end subroutine update_laser_particle_acc
#else
  end subroutine update_laser_particle
#endif
  
end module warpx_laser_module
