/* Copyright 2019 Axel Huebl, Maxence Thevenet, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef INJECTOR_DENSITY_H_
#define INJECTOR_DENSITY_H_

#include "Parser/GpuParser.H"
#include "CustomDensityProb.H"
#include "Utils/WarpXConst.H"

#include <AMReX_Gpu.H>
#include <AMReX_Dim3.H>

// struct whose getDensity returns constant density.
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

// struct whose getDensity returns local density computed from parser.
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

    // InjectorDensityParser constructs this GpuParser from WarpXParser.
    GpuParser<3> m_parser;
};

// struct whose getDensity returns local density computed from predefined profile.
struct InjectorDensityPredefined
{
    InjectorDensityPredefined (std::string const& a_species_name) noexcept;

    void clear ();

    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getDensity (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        // Choices for profile are:
        // - parabolic_channel
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

            // Longitudinal profile, normalized to 1
            if        ((z-z_start)>=0               and
                       (z-z_start)<ramp_up ) {
                n = 0.5*(1.-std::cos(MathConst::pi*(z-z_start)/ramp_up));
            } else if ((z-z_start)>=ramp_up         and
                       (z-z_start)< ramp_up+plateau ) {
                n = 1.;
            } else if ((z-z_start)>=ramp_up+plateau and
                       (z-z_start)< ramp_up+plateau+ramp_down) {
                n = 0.5*(1.+std::cos(MathConst::pi*((z-z_start)-ramp_up-plateau)/ramp_down));
            } else {
                n = 0.;
            }
            // Multiply by transverse profile, and physical density
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

// Base struct for density injector.
// InjectorDensity contains a union (called Object) that holds any one
// instance of:
// - InjectorDensityConstant  : to generate constant density;
// - InjectorDensityParser    : to generate density from parser;
// - InjectorDensityCustom    : to generate density from custom profile;
// - InjectorDensityPredefined: to generate density from predefined profile;
// The choice is made at runtime, depending in the constructor called.
// This mimics virtual functions, except the struct is stored in managed memory
// and member functions are made __host__ __device__ to run on CPU and GPU.
// This struct inherits from amrex::Gpu::Managed to provide new and delete
// operators in managed memory when running on GPU. Nothing special on CPU.
struct InjectorDensity
    : public amrex::Gpu::Managed
{
    // This constructor stores a InjectorDensityConstant in union object.
    InjectorDensity (InjectorDensityConstant* t, amrex::Real a_rho)
        : type(Type::constant),
          object(t,a_rho)
    { }

    // This constructor stores a InjectorDensityParser in union object.
    InjectorDensity (InjectorDensityParser* t, WarpXParser const& a_parser)
        : type(Type::parser),
          object(t,a_parser)
    { }

    // This constructor stores a InjectorDensityCustom in union object.
    InjectorDensity (InjectorDensityCustom* t, std::string const& a_species_name)
        : type(Type::custom),
          object(t,a_species_name)
    { }

    // This constructor stores a InjectorDensityPredefined in union object.
    InjectorDensity (InjectorDensityPredefined* t, std::string const& a_species_name)
        : type(Type::predefined),
          object(t,a_species_name)
    { }

    // Explicitly prevent the compiler from generating copy constructors
    // and copy assignment operators.
    InjectorDensity (InjectorDensity const&) = delete;
    InjectorDensity (InjectorDensity&&) = delete;
    void operator= (InjectorDensity const&) = delete;
    void operator= (InjectorDensity &&) = delete;

    ~InjectorDensity ();

    // call getDensity from the object stored in the union
    // (the union is called Object, and the instance is called object).
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

    // An instance of union Object constructs and stores any one of
    // the objects declared (constant or parser or custom or predefined).
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
