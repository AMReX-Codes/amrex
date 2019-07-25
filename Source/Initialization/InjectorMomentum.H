#ifndef INJECTOR_MOMENTUM_H_
#define INJECTOR_MOMENTUM_H_

#include <AMReX_Gpu.H>
#include <AMReX_Dim3.H>
#include <GpuParser.H>
#include <CustomMomentumProb.H>

struct InjectorMomentumConstant
{
    InjectorMomentumConstant (amrex::Real a_ux, amrex::Real a_uy, amrex::Real a_uz) noexcept
        : m_ux(a_ux), m_uy(a_uy), m_uz(a_uz) {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real, amrex::Real, amrex::Real) const noexcept
    {
        return amrex::XDim3{m_ux,m_uy,m_uz};
    }
private:
    amrex::Real m_ux, m_uy, m_uz;
};

struct InjectorMomentumGaussian
{
    InjectorMomentumGaussian (amrex::Real a_ux_m, amrex::Real a_uy_m,
                              amrex::Real a_uz_m, amrex::Real a_ux_th,
                              amrex::Real a_uy_th, amrex::Real a_uz_th) noexcept
        : m_ux_m(a_ux_m), m_uy_m(a_uy_m), m_uz_m(a_uz_m),
          m_ux_th(a_ux_th), m_uy_th(a_uy_th), m_uz_th(a_uz_th)
        {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return amrex::XDim3{amrex::RandomNormal(m_ux_m, m_ux_th),
                            amrex::RandomNormal(m_uy_m, m_uy_th),
                            amrex::RandomNormal(m_uz_m, m_uz_th)};
    }
private:
    amrex::Real m_ux_m, m_uy_m, m_uz_m;
    amrex::Real m_ux_th, m_uy_th, m_uz_th;
};

struct InjectorMomentumRadialExpansion
{
    InjectorMomentumRadialExpansion (amrex::Real a_u_over_r) noexcept
        : u_over_r(a_u_over_r)
        {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return {x*u_over_r, y*u_over_r, z*u_over_r};
    }

private:
    amrex::Real u_over_r;
};

struct InjectorMomentumParser
{
    InjectorMomentumParser (WarpXParser const& a_ux_parser,
                            WarpXParser const& a_uy_parser,
                            WarpXParser const& a_uz_parser) noexcept
        : m_ux_parser(a_ux_parser), m_uy_parser(a_uy_parser),
          m_uz_parser(a_uz_parser) {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return amrex::XDim3{m_ux_parser(x,y,z),m_uy_parser(x,y,z),m_uz_parser(x,y,z)};
    }

    GpuParser m_ux_parser, m_uy_parser, m_uz_parser;
};

struct InjectorMomentum
    : public amrex::Gpu::Managed
{
    InjectorMomentum (InjectorMomentumConstant* t,
                      amrex::Real a_ux, amrex::Real a_uy, amrex::Real a_uz)
        : type(Type::constant),
          object(t, a_ux, a_uy, a_uz)
    { }

    InjectorMomentum (InjectorMomentumParser* t,
                      WarpXParser const& a_ux_parser,
                      WarpXParser const& a_uy_parser,
                      WarpXParser const& a_uz_parser)
        : type(Type::parser),
          object(t, a_ux_parser, a_uy_parser, a_uz_parser)
    { }

    InjectorMomentum (InjectorMomentumGaussian* t,
                      amrex::Real a_ux_m, amrex::Real a_uy_m, amrex::Real a_uz_m,
                      amrex::Real a_ux_th, amrex::Real a_uy_th, amrex::Real a_uz_th)
        : type(Type::gaussian),
          object(t,a_ux_m,a_uy_m,a_uz_m,a_ux_th,a_uy_th,a_uz_th)
    { }

    InjectorMomentum (InjectorMomentumCustom* t,
                      std::string const& a_species_name)
        : type(Type::custom),
          object(t, a_species_name)
    { }

    InjectorMomentum (InjectorMomentumRadialExpansion* t,
                      amrex::Real u_over_r)
        : type(Type::radial_expansion),
          object(t, u_over_r)
    { }

    InjectorMomentum (InjectorMomentum const&) = delete;
    InjectorMomentum (InjectorMomentum&&) = delete;
    void operator= (InjectorMomentum const&) = delete;
    void operator= (InjectorMomentum &&) = delete;

    ~InjectorMomentum ();

    std::size_t sharedMemoryNeeded () const noexcept;

    bool useRandom () const noexcept;

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        switch (type)
        {
        case Type::parser:
        {
            return object.parser.getMomentum(x,y,z);
        }
        case Type::gaussian:
        {
            return object.gaussian.getMomentum(x,y,z);
        }
        case Type::constant:
        {
            return object.constant.getMomentum(x,y,z);
        }
        case Type::radial_expansion:
        {
            return object.radial_expansion.getMomentum(x,y,z);
        }
        case Type::custom:
        {
            return object.custom.getMomentum(x,y,z);
        }
        default:
        {
            amrex::Abort("InjectorMomentum: unknown type");
            return {0.0,0.0,0.0};
        }
        }
    }

private:
    enum struct Type { constant, custom, gaussian, radial_expansion, parser };
    Type type;

    union Object {
        Object (InjectorMomentumConstant*,
                amrex::Real a_ux, amrex::Real a_uy, amrex::Real a_uz) noexcept
            : constant(a_ux,a_uy,a_uz) {}
        Object (InjectorMomentumCustom*,
                std::string const& a_species_name) noexcept
            : custom(a_species_name) {}
        Object (InjectorMomentumGaussian*,
                amrex::Real a_ux_m, amrex::Real a_uy_m,
                amrex::Real a_uz_m, amrex::Real a_ux_th,
                amrex::Real a_uy_th, amrex::Real a_uz_th) noexcept
            : gaussian(a_ux_m,a_uy_m,a_uz_m,a_ux_th,a_uy_th,a_uz_th) {}
        Object (InjectorMomentumRadialExpansion*,
                amrex::Real u_over_r) noexcept
            : radial_expansion(u_over_r) {}
        Object (InjectorMomentumParser*,
                WarpXParser const& a_ux_parser,
                WarpXParser const& a_uy_parser,
                WarpXParser const& a_uz_parser) noexcept
            : parser(a_ux_parser, a_uy_parser, a_uz_parser) {}
        InjectorMomentumConstant constant;
        InjectorMomentumCustom   custom;
        InjectorMomentumGaussian gaussian;
        InjectorMomentumRadialExpansion radial_expansion;
        InjectorMomentumParser   parser;
    };
    Object object;
};

#endif
