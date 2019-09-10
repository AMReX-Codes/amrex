
#include <WarpX_Complex.H>
#include <LaserParticleContainer.H>

using namespace amrex;

/* \brief compute field amplitude for a Gaussian laser, at particles' position
 *
 * Both Xp and Yp are given in laser plane coordinate.
 * For each particle with position Xp and Yp, this routine computes the
 * amplitude of the laser electric field, stored in array amplitude.
 * 
 * \param np: number of laser particles
 * \param Xp: pointer to first component of positions of laser particles
 * \param Yp: pointer to second component of positions of laser particles
 * \param t: Current physical time
 * \param amplitude: pointer to array of field amplitude.
 */
void 
LaserParticleContainer::gaussian_laser_profile (
    const int np, Real const * const Xp, Real const * const Yp,
    Real t, Real * const amplitude)
{
    Complex I(0,1);
    // Calculate a few factors which are independent of the macroparticle
    const Real k0 = 2.*MathConst::pi/wavelength;
    const Real inv_tau2 = 1. / (profile_duration * profile_duration);
    const Real oscillation_phase = k0 * PhysConst::c * ( t - profile_t_peak );
    // The coefficients below contain info about Gouy phase,
    // laser diffraction, and phase front curvature
    const Complex diffract_factor = 1. + I * profile_focal_distance
        * 2./( k0 * profile_waist * profile_waist );
    const Complex inv_complex_waist_2 = 1./( profile_waist*profile_waist * diffract_factor );

    // Time stretching due to STCs and phi2 complex envelope
    // (1 if zeta=0, beta=0, phi2=0)
    const Complex stretch_factor = 1. + 4. * 
        (zeta+beta*profile_focal_distance) * (zeta+beta*profile_focal_distance)
        * (inv_tau2*inv_complex_waist_2) + 
        2.*I*(phi2 - beta*beta*k0*profile_focal_distance) * inv_tau2;

    // Amplitude and monochromatic oscillations
    Complex prefactor = e_max * MathFunc::exp( I * oscillation_phase );

    // Because diffract_factor is a complex, the code below takes into
    // account the impact of the dimensionality on both the Gouy phase
    // and the amplitude of the laser
#if (AMREX_SPACEDIM == 3)
    prefactor = prefactor / diffract_factor;
#elif (AMREX_SPACEDIM == 2)
    prefactor = prefactor / MathFunc::sqrt(diffract_factor);
#endif

    // Copy member variables to tmp copies for GPU runs.
    Real tmp_profile_t_peak = profile_t_peak;
    Real tmp_beta = beta;
    Real tmp_zeta = zeta;
    // Loop through the macroparticle to calculate the proper amplitude
    amrex::ParallelFor(
        np, 
        [=] AMREX_GPU_DEVICE (int i) {
            const Complex stc_exponent = 1./stretch_factor * inv_tau2 *
                MathFunc::pow((t - tmp_profile_t_peak - 
                               tmp_beta*k0*(Xp[i]*std::cos(theta_stc) + Yp[i]*std::sin(theta_stc)) -
                               2.*I*(Xp[i]*std::cos(theta_stc) + Yp[i]*std::sin(theta_stc))
                               *( tmp_zeta - tmp_beta*profile_focal_distance ) * inv_complex_waist_2),2);
            // stcfactor = everything but complex transverse envelope
            const Complex stcfactor = prefactor * MathFunc::exp( - stc_exponent );
            // Exp argument for transverse envelope
            const Complex exp_argument = - ( Xp[i]*Xp[i] + Yp[i]*Yp[i] ) * inv_complex_waist_2;
            // stcfactor + transverse envelope
            amplitude[i] = ( stcfactor * MathFunc::exp( exp_argument ) ).real();
        }
        );
}

/* \brief compute field amplitude for a Harris laser function, at particles' position
 *
 * Both Xp and Yp are given in laser plane coordinate.
 * For each particle with position Xp and Yp, this routine computes the
 * amplitude of the laser electric field, stored in array amplitude.
 * 
 * \param np: number of laser particles
 * \param Xp: pointer to first component of positions of laser particles
 * \param Yp: pointer to second component of positions of laser particles
 * \param t: Current physical time
 * \param amplitude: pointer to array of field amplitude.
 */
void 
LaserParticleContainer::harris_laser_profile (
    const int np, Real const * const Xp, Real const * const Yp,
    Real t, Real * const amplitude)
{
    // This function uses the Harris function as the temporal profile of the pulse
    const Real omega0 = 2.*MathConst::pi*PhysConst::c/wavelength;
    const Real zR = MathConst::pi * profile_waist*profile_waist / wavelength;
    const Real wz = profile_waist * 
        std::sqrt(1. + profile_focal_distance*profile_focal_distance/zR*zR);
    const Real inv_wz_2 = 1./(wz*wz);
    Real inv_Rz;
    if (profile_focal_distance == 0.){ 
        inv_Rz = 0.;
    } else {
        inv_Rz = -profile_focal_distance / 
            ( profile_focal_distance*profile_focal_distance + zR*zR );
    }
    const Real arg_env = 2.*MathConst::pi*t/profile_duration;

    // time envelope is given by the Harris function
    Real time_envelope = 0.;
    if (t < profile_duration)
        time_envelope = 1./32. * (10. - 15.*std::cos(arg_env) + 
                                  6.*std::cos(2.*arg_env) - 
                                  std::cos(3.*arg_env));

    // Copy member variables to tmp copies for GPU runs.
    Real tmp_e_max = e_max;
    // Loop through the macroparticle to calculate the proper amplitude
    amrex::ParallelFor(
        np,
        [=] AMREX_GPU_DEVICE (int i) {
            const Real space_envelope = 
                std::exp(- ( Xp[i]*Xp[i] + Yp[i]*Yp[i] ) * inv_wz_2);
            const Real arg_osc = omega0*t - omega0/PhysConst::c*
                (Xp[i]*Xp[i] + Yp[i]*Yp[i]) * inv_Rz / 2.;
            const Real oscillations = std::cos(arg_osc);
            amplitude[i] = tmp_e_max * time_envelope *
                space_envelope * oscillations;
        }
        );
}
