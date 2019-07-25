#ifndef INJECTOR_DENSITY_H_
#define INJECTOR_DENSITY_H_

#include <AMReX_Gpu.H>
#include <AMReX_Dim3.H>
#include <GpuParser.H>
#include <CustomDensityProb.H>
#include <WarpXConst.H>

struct InjectorDensityConstant
{
    InjectorDensityConstant (amrex::Real a_rho) noexcept : m_rho(a_rho) {}

    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getDensity (amrex::Real, amrex::Real, amrex::Real) const noexcept
    {
        return m_rho;
    }

private:
    amrex::Real m_rho;
};

struct InjectorDensityParser
{
    InjectorDensityParser (WarpXParser const& a_parser) noexcept
        : m_parser(a_parser) {}

    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getDensity (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return m_parser(x,y,z);
    }

    GpuParser m_parser;
};

struct InjectorDensityPredefined
{
    InjectorDensityPredefined (std::string const& a_species_name) noexcept;

    void clear ();

    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getDensity (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        switch (profile)
        {
        case Profile::parabolic_channel:
        {
            amrex::Real z_start   = p[0];
            amrex::Real ramp_up   = p[1];
            amrex::Real plateau   = p[2];
            amrex::Real ramp_down = p[3];
            amrex::Real rc        = p[4];
            amrex::Real n0        = p[5];
            amrex::Real n;
            amrex::Real kp = PhysConst::q_e/PhysConst::c
                *std::sqrt( n0/(PhysConst::m_e*PhysConst::ep0) );

            if        ((z-z_start)>=0               and
                       (z-z_start)<ramp_up ) {
                n = (z-z_start)/ramp_up;
            } else if ((z-z_start)>=ramp_up         and
                       (z-z_start)< ramp_up+plateau ) {
                n = 1.;
            } else if ((z-z_start)>=ramp_up+plateau and
                       (z-z_start)< ramp_up+plateau+ramp_down) {
                n = 1.-((z-z_start)-ramp_up-plateau)/ramp_down;
            } else {
                n = 0.;
            }
            n *= n0*(1.+4.*(x*x+y*y)/(kp*kp*rc*rc*rc*rc));
            return n;
        }
        default:
            amrex::Abort("InjectorDensityPredefined: how did we get here?");
            return 0.0;
        }
    }

private:
    enum struct Profile { null, parabolic_channel };
    Profile profile;
    amrex::Real* p;
};

struct InjectorDensity
    : public amrex::Gpu::Managed
{
    InjectorDensity (InjectorDensityConstant* t, amrex::Real a_rho)
        : type(Type::constant),
          object(t,a_rho)
    { }

    InjectorDensity (InjectorDensityParser* t, WarpXParser const& a_parser)
        : type(Type::parser),
          object(t,a_parser)
    { }

    InjectorDensity (InjectorDensityCustom* t, std::string const& a_species_name)
        : type(Type::custom),
          object(t,a_species_name)
    { }

    InjectorDensity (InjectorDensityPredefined* t, std::string const& a_species_name)
        : type(Type::predefined),
          object(t,a_species_name)
    { }

    InjectorDensity (InjectorDensity const&) = delete;
    InjectorDensity (InjectorDensity&&) = delete;
    void operator= (InjectorDensity const&) = delete;
    void operator= (InjectorDensity &&) = delete;

    ~InjectorDensity ();

    std::size_t sharedMemoryNeeded () const noexcept;
    bool useRandom () const noexcept { return false; }

    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getDensity (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        switch (type)
        {
        case Type::parser:
        {
            return object.parser.getDensity(x,y,z);
        }
        case Type::constant:
        {
            return object.constant.getDensity(x,y,z);
        }
        case Type::custom:
        {
            return object.custom.getDensity(x,y,z);
        }
        case Type::predefined:
        {
            return object.predefined.getDensity(x,y,z);
        }
        default:
        {
            amrex::Abort("InjectorDensity: unknown type");
            return 0.0;
        }
        }
    }

private:
    enum struct Type { constant, custom, predefined, parser };
    Type type;

    union Object {
        Object (InjectorDensityConstant*, amrex::Real a_rho) noexcept
            : constant(a_rho) {}
        Object (InjectorDensityParser*, WarpXParser const& a_parser) noexcept
            : parser(a_parser) {}
        Object (InjectorDensityCustom*, std::string const& a_species_name) noexcept
            : custom(a_species_name) {}
        Object (InjectorDensityPredefined*, std::string const& a_species_name) noexcept
            : predefined(a_species_name) {}
        InjectorDensityConstant   constant;
        InjectorDensityParser     parser;
        InjectorDensityCustom     custom;
        InjectorDensityPredefined predefined;
    };
    Object object;
};

#endif
