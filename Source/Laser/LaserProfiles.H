#ifndef WARPX_LaserProfiles_H_
#define WARPX_LaserProfiles_H_

#include <AMReX_REAL.H>
#include <AMReX_ParamParse.H>
#include <AMReX_Vector.H>

#include <WarpXParser.H>

#include <map>
#include <string>
#include <memory>
#include <functional>
#include <limits>

//_________________________________________Abstract interface for laser profile classes

class ILaserProfile
{
public:
    virtual void
    init (const amrex::ParmParse& pp) = 0;

    virtual void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) = 0;

    virtual ~ILaserProfile(){};
};

//_________________________________________Declaration of laser profile classes

class GaussianLaserProfile : public ILaserProfile
{

public:
    void init (const amrex::ParmParse& pp) override final;

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
        amrex::Vector<amrex::Real> stc_direction; //What should I do?
    } m_params;
};

class HarrisLaserProfile : public ILaserProfile
{

public:
    void init (const amrex::ParmParse& pp) override final;

    void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) override final;

private:
    struct {
        amrex::Real waist;
        amrex::Real duration;
        amrex::Real focal_distance;
    } m_params;
};

class FieldFunctionLaserProfile : public ILaserProfile
{

public:
    void init (const amrex::ParmParse& pp) override final;

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

//_________________________________________Map of laser profile names --> laser profile classes

std::map<
std::string,
std::function<std::unique_ptr<ILaserProfile>()>
>
laser_profiles_dictionary =
{
    {"gaussian", [] () {return std::make_unique<GaussianLaserProfile>();} },
    {"harris", [] () {return std::make_unique<HarrisLaserProfile>();} },
    {"parse_field_function", [] () {return std::make_unique<FieldFunctionLaserProfile>();} }
};

#endif //WARPX_LaserProfiles_H_
