/* Copyright 2019 Axel Huebl, Cameron Yang, Maxence Thevenet
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef INJECTOR_MOMENTUM_H_
#define INJECTOR_MOMENTUM_H_

#include "CustomMomentumProb.H"
#include "Parser/GpuParser.H"

#include <AMReX_Gpu.H>
#include <AMReX_Dim3.H>

// struct whose getMomentum returns constant momentum.
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

// struct whose getMomentum returns momentum for 1 particle, from random
// gaussian distribution.
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

// struct whose getMomentum returns momentum for 1 particle with relativistic
// drift velocity beta, from the Maxwell-Boltzmann distribution.
struct InjectorMomentumBoltzmann
{
    // Constructor whose inputs are:
    // the temperature parameter theta,
    // boost velocity/c beta,
    // and boost direction dir respectively.
    InjectorMomentumBoltzmann(amrex::Real t, amrex::Real b, int d) noexcept
        : dir(d), beta(b), vave(std::sqrt(2.*t))
        {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        amrex::Real x1, x2, gamma;
        amrex::Real u[3];
        x1 = amrex::Random();
        x2 = amrex::Random();
        // Each value of sqrt(-log(x1))*sin(2*pi*x2) is a sample from a Gaussian
        // distribution with sigma = average velocity / c
        // using the Box-Mueller Method.
        u[(dir+1)%3] = vave*std::sqrt(-std::log(x1)) *std::sin(2*M_PI*x2);
        u[(dir+2)%3] = vave*std::sqrt(-std::log(x1)) *std::cos(2*M_PI*x2);
        u[dir] = vave*std::sqrt(-std::log(amrex::Random()))*
            std::sin(2*M_PI*amrex::Random());
        gamma = std::pow(u[0],2)+std::pow(u[1],2)+std::pow(u[2],2);
        gamma = std::sqrt(1+gamma);
        // The following condition is equtaion 32 in Zenitani 2015
        // (Phys. Plasmas 22, 042116) , called the flipping method. It
        // transforms the intergral: d3x' -> d3x  where d3x' is the volume
        // element for positions in the boosted frame. The particle positions
        // and densities can be initialized in the simulation frame.
        // The flipping method can transform any symmetric distribution from one
        // reference frame to another moving at a relative velocity of beta.
        // An equivalent alternative to this method native to WarpX would be to
        // initialize the particle positions and densities in the frame moving
        // at speed beta, and then perform a Lorentz transform on the positions
        // and MB sampled velocities to the simulation frame.
        x1 = amrex::Random();
        if(-beta*u[dir]/gamma > x1)
        {
          u[dir] = -u[dir];
        }
        // This Lorentz transform is equation 17 in Zenitani.
        // It transforms the integral d3u' -> d3u
        // where d3u' is the volume element for momentum in the boosted frame.
        u[dir] = 1/std::sqrt(1-pow(beta,2))*(u[dir]+gamma*beta);
        // Note that if beta = 0 then the flipping method and Lorentz transform
        // have no effect on the u[dir] direction.
        return amrex::XDim3 {u[0],u[1],u[2]};
    }

private:
    int dir;
    amrex::Real beta, vave;
};

// struct whose getMomentum returns momentum for 1 particle with relativistc
// drift velocity beta, from the Maxwell-Juttner distribution. Method is from
// Zenitani 2015 (Phys. Plasmas 22, 042116).
struct InjectorMomentumJuttner
{
    // Constructor whose inputs are:
    // the temperature parameter theta,
    // boost velocity/c beta,
    // and boost direction dir respectively.
    InjectorMomentumJuttner(amrex::Real t, amrex::Real b, int d) noexcept
        : dir(d), beta(b), theta(t)
        {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        // Sobol method for sampling MJ Speeds,
        // from Zenitani 2015 (Phys. Plasmas 22, 042116).
        amrex::Real x1, x2, gamma;
        amrex::Real u [3];
        x1 = 0.;
        gamma = 0.;
        u[dir] = 0.;
        // This condition is equation 10 in Zenitani,
        // though x1 is defined differently.
        while(u[dir]-gamma <= x1)
        {
            u[dir] = -theta*
                std::log(amrex::Random()*amrex::Random()*amrex::Random());
            gamma = std::sqrt(1+std::pow(u[dir],2));
            x1 = theta*std::log(amrex::Random());
        }
        // The following code samples a random unit vector
        // and multiplies the result by speed u[dir].
        x1 = amrex::Random();
        x2 = amrex::Random();
        // Direction dir is an input parameter that sets the boost direction:
        // 'x' -> d = 0, 'y' -> d = 1, 'z' -> d = 2.
        u[(dir+1)%3] = 2*u[dir]*std::sqrt(x1*(1-x1))*std::sin(2*M_PI*x2);
        u[(dir+2)%3] = 2*u[dir]*std::sqrt(x1*(1-x1))*std::cos(2*M_PI*x2);
        // The value of dir is the boost direction to be transformed.
        u[dir] = u[dir]*(2*x1-1);
        x1 = amrex::Random();
        // The following condition is equtaion 32 in Zenitani, called
        // The flipping method. It transforms the intergral: d3x' -> d3x
        // where d3x' is the volume element for positions in the boosted frame.
        // The particle positions and densities can be initialized in the
        // simulation frame with this method.
        // The flipping method can similarly transform any
        // symmetric distribution from one reference frame to another moving at
        // a relative velocity of beta.
        // An equivalent alternative to this method native to WarpX
        // would be to initialize the particle positions and densities in the
        // frame moving at speed beta, and then perform a Lorentz transform
        // on their positions and MJ sampled velocities to the simulation frame.
        if(-beta*u[dir]/gamma>x1)
        {
            u[dir] = -u[dir];
        }
        // This Lorentz transform is equation 17 in Zenitani.
        // It transforms the integral d3u' -> d3u
        // where d3u' is the volume element for momentum in the boosted frame.
        u[dir] = 1/std::sqrt(1-pow(beta,2))*(u[dir]+gamma*beta);
        // Note that if beta = 0 then the flipping method and Lorentz transform
        // have no effect on the u[dir] direction.
        return amrex::XDim3 {u[0],u[1],u[2]};
    }

private:
    int dir;
    amrex::Real beta, theta;
};


// struct whose getMomentum returns momentum for 1 particle, for
// radial expansion
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

// struct whose getMomentumm returns local momentum computed from parser.
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

    GpuParser<3> m_ux_parser, m_uy_parser, m_uz_parser;
};

// Base struct for momentum injector.
// InjectorMomentum contains a union (called Object) that holds any one
// instance of:
// - InjectorMomentumConstant       : to generate constant density;
// - InjectorMomentumGaussian       : to generate gaussian distribution;
// - InjectorMomentumRadialExpansion: to generate radial expansion;
// - InjectorMomentumParser         : to generate momentum from parser;
// The choice is made at runtime, depending in the constructor called.
// This mimics virtual functions, except the struct is stored in managed memory
// and member functions are made __host__ __device__ to run on CPU and GPU.
// This struct inherits from amrex::Gpu::Managed to provide new and delete
// operators in managed memory when running on GPU. Nothing special on CPU.
struct InjectorMomentum
    : public amrex::Gpu::Managed
{
    // This constructor stores a InjectorMomentumConstant in union object.
    InjectorMomentum (InjectorMomentumConstant* t,
                      amrex::Real a_ux, amrex::Real a_uy, amrex::Real a_uz)
        : type(Type::constant),
          object(t, a_ux, a_uy, a_uz)
    { }

    // This constructor stores a InjectorMomentumParser in union object.
    InjectorMomentum (InjectorMomentumParser* t,
                      WarpXParser const& a_ux_parser,
                      WarpXParser const& a_uy_parser,
                      WarpXParser const& a_uz_parser)
        : type(Type::parser),
          object(t, a_ux_parser, a_uy_parser, a_uz_parser)
    { }

    // This constructor stores a InjectorMomentumGaussian in union object.
    InjectorMomentum (InjectorMomentumGaussian* t,
                      amrex::Real a_ux_m, amrex::Real a_uy_m, amrex::Real a_uz_m,
                      amrex::Real a_ux_th, amrex::Real a_uy_th, amrex::Real a_uz_th)
        : type(Type::gaussian),
          object(t,a_ux_m,a_uy_m,a_uz_m,a_ux_th,a_uy_th,a_uz_th)
    { }

    InjectorMomentum (InjectorMomentumBoltzmann* t,
                       amrex::Real theta, amrex::Real beta, int dir)
         : type(Type::boltzmann),
           object(t, theta, beta, dir)
    { }

     // This constructor stores a InjectorMomentumJuttner in union object.
     InjectorMomentum (InjectorMomentumJuttner* t,
                       amrex::Real theta, amrex::Real beta, int dir)
         : type(Type::juttner),
           object(t, theta, beta, dir)
    { }

    // This constructor stores a InjectorMomentumCustom in union object.
    InjectorMomentum (InjectorMomentumCustom* t,
                      std::string const& a_species_name)
        : type(Type::custom),
          object(t, a_species_name)
    { }

    // This constructor stores a InjectorMomentumRadialExpansion in union object.
    InjectorMomentum (InjectorMomentumRadialExpansion* t,
                      amrex::Real u_over_r)
        : type(Type::radial_expansion),
          object(t, u_over_r)
    { }

    // Explicitly prevent the compiler from generating copy constructors
    // and copy assignment operators.
    InjectorMomentum (InjectorMomentum const&) = delete;
    InjectorMomentum (InjectorMomentum&&) = delete;
    void operator= (InjectorMomentum const&) = delete;
    void operator= (InjectorMomentum &&) = delete;

    ~InjectorMomentum ();

    // call getMomentum from the object stored in the union
    // (the union is called Object, and the instance is called object).
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
        case Type::boltzmann:
        {
            return object.boltzmann.getMomentum(x,y,z);
        }
        case Type::juttner:
        {
            return object.juttner.getMomentum(x,y,z);
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
    enum struct Type { constant, custom, gaussian, boltzmann, juttner, radial_expansion, parser};
    Type type;

    // An instance of union Object constructs and stores any one of
    // the objects declared (constant or custom or gaussian or
    // radial_expansion or parser).
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
        Object (InjectorMomentumBoltzmann*,
                amrex::Real t, amrex::Real b, int dir) noexcept
            : boltzmann(t,b,dir) {}
        Object (InjectorMomentumJuttner*,
                amrex::Real t, amrex::Real b, int dir) noexcept
            : juttner(t,b,dir) {}
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
        InjectorMomentumBoltzmann boltzmann;
        InjectorMomentumJuttner juttner;
        InjectorMomentumRadialExpansion radial_expansion;
        InjectorMomentumParser   parser;
    };
    Object object;
};

#endif
