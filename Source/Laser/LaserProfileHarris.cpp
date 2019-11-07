#include <LaserProfiles.H>

#include <WarpX_Complex.H>

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
HarrisLaserProfile::fill_amplitude (
    const int np, Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude)
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
