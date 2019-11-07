#include <LaserProfiles.H>

#include <WarpX_Complex.H>

using namespace amrex;

void
GaussianLaserProfile::init (const amrex::ParmParse& pp)
{
    // Parse the properties of the Gaussian profile
    pp.get("profile_waist", m_params.profile_waist);
    pp.get("profile_duration", m_params.profile_duration);
    pp.get("profile_t_peak", m_params.profile_t_peak);
    pp.get("profile_focal_distance", m_params.profile_focal_distance);
    //stc_direction = p_X; //TO HANDLE!
    pp.queryarr("stc_direction", m_params.stc_direction);
    pp.query("zeta", m_params.zeta);
    pp.query("beta", m_params.beta);
    pp.query("phi2", m_params.phi2);
}

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
GaussianLaserProfile::fill_amplitude (
    const int np, Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude)
{
    Complex I(0,1);
    // Calculate a few factors which are independent of the macroparticle
    const Real k0 = 2.*MathConst::pi/wavelength;
    const Real inv_tau2 = 1. / (profile_duration * profile_duration);
    const Real oscillation_phase = k0 * PhysConst::c * ( t - profile_t_peak );
    // The coefficients below contain info about Gouy phase,
    // laser diffraction, and phase front curvature
    const Complex diffract_factor = 1._rt + I * profile_focal_distance
        * 2._rt/( k0 * profile_waist * profile_waist );
    const Complex inv_complex_waist_2 = 1._rt / ( profile_waist*profile_waist * diffract_factor );

    // Time stretching due to STCs and phi2 complex envelope
    // (1 if zeta=0, beta=0, phi2=0)
    const Complex stretch_factor = 1._rt + 4._rt *
        (zeta+beta*profile_focal_distance) * (zeta+beta*profile_focal_distance)
        * (inv_tau2*inv_complex_waist_2) +
        2._rt *I*(phi2 - beta*beta*k0*profile_focal_distance) * inv_tau2;

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
    Real tmp_theta_stc = theta_stc;
    Real tmp_profile_focal_distance = profile_focal_distance;
    // Loop through the macroparticle to calculate the proper amplitude
    amrex::ParallelFor(
        np,
        [=] AMREX_GPU_DEVICE (int i) {
            const Complex stc_exponent = 1._rt / stretch_factor * inv_tau2 *
                MathFunc::pow((t - tmp_profile_t_peak -
                               tmp_beta*k0*(Xp[i]*std::cos(tmp_theta_stc) + Yp[i]*std::sin(tmp_theta_stc)) -
                               2._rt *I*(Xp[i]*std::cos(tmp_theta_stc) + Yp[i]*std::sin(tmp_theta_stc))
                               *( tmp_zeta - tmp_beta*tmp_profile_focal_distance ) * inv_complex_waist_2),2);
            // stcfactor = everything but complex transverse envelope
            const Complex stcfactor = prefactor * MathFunc::exp( - stc_exponent );
            // Exp argument for transverse envelope
            const Complex exp_argument = - ( Xp[i]*Xp[i] + Yp[i]*Yp[i] ) * inv_complex_waist_2;
            // stcfactor + transverse envelope
            amplitude[i] = ( stcfactor * MathFunc::exp( exp_argument ) ).real();
        }
        );
}
