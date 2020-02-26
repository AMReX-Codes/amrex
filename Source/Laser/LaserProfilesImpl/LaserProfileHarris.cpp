/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Laser/LaserProfiles.H"
#include "Utils/WarpX_Complex.H"
#include "Utils/WarpXConst.H"


using namespace amrex;

void
WarpXLaserProfiles::HarrisLaserProfile::init (
    const amrex::ParmParse& ppl,
    const amrex::ParmParse& /* ppc */,
    CommonLaserParameters params)
{
    // Parse the properties of the Harris profile
    ppl.get("profile_waist", m_params.waist);
    ppl.get("profile_duration", m_params.duration);
    ppl.get("profile_focal_distance", m_params.focal_distance);
    //Copy common params
    m_common_params = params;
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
WarpXLaserProfiles::HarrisLaserProfile::fill_amplitude (
    const int np, Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude) const
{
    // This function uses the Harris function as the temporal profile of the pulse
    const Real omega0 =
        2._rt*MathConst::pi*PhysConst::c/m_common_params.wavelength;
    const Real zR = MathConst::pi * m_params.waist*m_params.waist
        / m_common_params.wavelength;
    const Real wz = m_params.waist *
        std::sqrt(1._rt + m_params.focal_distance*m_params.focal_distance/(zR*zR));
    const Real inv_wz_2 = 1._rt/(wz*wz);
    Real inv_Rz;
    if (m_params.focal_distance == 0.){
        inv_Rz = 0.;
    } else {
        inv_Rz = -m_params.focal_distance /
            ( m_params.focal_distance*m_params.focal_distance + zR*zR );
    }
    const Real arg_env = 2._rt*MathConst::pi*t/m_params.duration;

    // time envelope is given by the Harris function
    Real time_envelope = 0.;
    if (t < m_params.duration)
        time_envelope = 1._rt/32._rt * (10._rt - 15._rt*std::cos(arg_env) +
                                  6._rt*std::cos(2._rt*arg_env) -
                                  std::cos(3._rt*arg_env));

    // Copy member variables to tmp copies for GPU runs.
    const auto tmp_e_max = m_common_params.e_max;
    // Loop through the macroparticle to calculate the proper amplitude
    amrex::ParallelFor(
        np,
        [=] AMREX_GPU_DEVICE (int i) {
            const Real space_envelope =
                std::exp(- ( Xp[i]*Xp[i] + Yp[i]*Yp[i] ) * inv_wz_2);
            const Real arg_osc = omega0*t - omega0/PhysConst::c*
                (Xp[i]*Xp[i] + Yp[i]*Yp[i]) * inv_Rz / 2._rt;
            const Real oscillations = std::cos(arg_osc);
            amplitude[i] = tmp_e_max * time_envelope *
                space_envelope * oscillations;
        }
        );
}
