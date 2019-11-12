#ifndef WARPX_LaserProfiles_H_
#define WARPX_LaserProfiles_H_

#include <AMReX_REAL.H>
#include <WarpXParser.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <WarpXParser.H>

#include <map>
#include <string>
#include <memory>
#include <functional>
#include <limits>

namespace WarpXLaserProfiles {

/** Common laser profile parameters
 *
 * Parameters for each laser profile as shared among all laser profile classes.
 */
struct CommonLaserParameters
{
    amrex::Real wavelength; //! central wavelength
    amrex::Real e_max;  //! maximum electric field at peak
    amrex::Vector<amrex::Real> p_X;// ! Polarization
    amrex::Vector<amrex::Real> nvec; //! Normal of the plane of the antenna
};


/** Abstract interface for laser profile classes
 *
 * Each new laser profile should inherit this interface and implement two
 * methods: init and fill_amplitude (described below).
 *
 * The declaration of a LaserProfile class should be placed in this file,
 * while the implementation of the methods should be in a dedicated file in
 * LaserProfilesImpl folder. LaserProfile classes should appear in
 * laser_profiles_dictionary to be used by LaserParticleContainer.
 */
class ILaserProfile
{
public:
    /** Initialize Laser Profile
     *
     * Reads the section of the inputfile relative to the laser beam
     * (e.g. <laser_name>.profile_t_peak, <laser_name>.profile_duration...)
     * and the "my_constants" section. It also receives some common
     * laser profile parameters. It uses these data to initialize the
     * member variables of the laser profile class.
     *
     * @param[in] ppl should be amrex::ParmParse(laser_name)
     * @param[in] ppc should be amrex::ParmParse("my_constants")
     * @param[in] params common laser profile parameters
     */
    virtual void
    init (
        const amrex::ParmParse& ppl,
        const amrex::ParmParse& ppc,
        CommonLaserParameters params) = 0;

    /** Fill Electric Field Amplitude for each particle of the antenna.
     *
     * Xp, Yp and amplitude must be arrays with the same length
     *
     * @param[in] Xp X coordinate of the particles of the antenna
     * @param[in] Yp Y coordinate of the particles of the antenna
     * @param[in] t time (seconds)
     * @param[out] amplitude of the electric field (V/m)
     */
    virtual void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) = 0;

    virtual ~ILaserProfile(){};
};

/**
 * Gaussian laser profile
 */
class GaussianLaserProfile : public ILaserProfile
{

public:
    void
    init (
        const amrex::ParmParse& ppl,
        const amrex::ParmParse& ppc,
        CommonLaserParameters params) override final;

    void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) override final;

private:
    struct {
        amrex::Real waist          = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real duration       = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real t_peak         = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real focal_distance = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real zeta           = 0;
        amrex::Real beta           = 0;
        amrex::Real phi2           = 0;

        amrex::Vector<amrex::Real> stc_direction; //! Direction of the spatio-temporal couplings
        amrex::Real theta_stc; //! Angle between polarization (p_X) and direction of spatiotemporal coupling (stc_direction)
    } m_params;

    CommonLaserParameters m_common_params;
};

/**
 * Harris laser profile
 */
class HarrisLaserProfile : public ILaserProfile
{

public:
    void
    init (
        const amrex::ParmParse& ppl,
        const amrex::ParmParse& ppc,
        CommonLaserParameters params) override final;

    void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) override final;

private:
    struct {
        amrex::Real waist          = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real duration       = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real focal_distance = std::numeric_limits<amrex::Real>::quiet_NaN();
    } m_params;

    CommonLaserParameters m_common_params;
};

/**
 * Laser profile defined by the used with an analytical expression
 */
class FieldFunctionLaserProfile : public ILaserProfile
{

public:
    void
    init (
        const amrex::ParmParse& ppl,
        const amrex::ParmParse& ppc,
        CommonLaserParameters params) override final;

    void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) override final;

private:
    struct{
        std::string field_function;
    } m_params;

     WarpXParser m_parser;
};

/**
 * Maps laser profile names to lambdas returing unique pointers
 * to the corresponding laser profile objects.
 */
const
std::map<
std::string,
std::function<std::unique_ptr<ILaserProfile>()>
>
laser_profiles_dictionary =
{
    {"gaussian",
        [] () {return std::make_unique<GaussianLaserProfile>();} },
    {"harris",
        [] () {return std::make_unique<HarrisLaserProfile>();} },
    {"parse_field_function",
        [] () {return std::make_unique<FieldFunctionLaserProfile>();} }
};

} //WarpXLaserProfiles

#endif //WARPX_LaserProfiles_H_
