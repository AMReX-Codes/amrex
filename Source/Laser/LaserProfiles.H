#ifndef WARPX_LaserProfiles_H_
#define WARPX_LaserProfiles_H_

#include <AMReX_REAL.H>
#include <AMReX_ParamParse.H>

#include <WarpXParser.H>

#include <map>
#include <string>
#include <memory>
#include <functional>

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
    struct params{
        amrex::Real waist;
        amrex::Real duration;
        amrex::Real t_peak;
        amrex::Real focal_distance;
        amrex::Real stc_direction;
        amrex::Real zeta;
        amrex::Real beta;
        amrex::Real phi2;
    };
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
    struct params{
        amrex::Real waist;
        amrex::Real duration;
        amrex::Real focal_distance;
    };
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
    struct params{
        std::string field_function;
    };

     WarpXParser parser;
};

//_________________________________________Map of laser profile names --> laser profile classes

const std::map<
std::string,
std::function<std::unique_ptr<ILaserProfile>>(const amrex::ParmParse&)
>
laser_profiles_dictionary =
{
    {"gaussian", [] (const amrex::ParmParse& pp) {
        return std::make_unique<GaussianLaserProfile>(pp)} },
    {"harris", [] (const amrex::ParmParse& pp) {
        return std::make_unique<HarrisLaserProfile>(pp)} },
    {"parse_field_function", [] (const amrex::ParmParse& pp) {
        return std::make_unique<FieldFunctionLaserProfile>(pp)} }
};

#endif //WARPX_LaserProfiles_H_
