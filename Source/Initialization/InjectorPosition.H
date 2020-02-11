/* Copyright 2019 Axel Huebl, David Grote, Maxence Thevenet
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef INJECTOR_POSITION_H_
#define INJECTOR_POSITION_H_

#include <AMReX_Gpu.H>
#include <AMReX_Dim3.H>
#include <AMReX_Utility.H>

// struct whose getPositionUnitBox returns x, y and z for a particle with
// random distribution inside a unit cell.
struct InjectorPositionRandom
{
    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getPositionUnitBox (int i_part, int ref_fac=1) const noexcept
    {
        return amrex::XDim3{amrex::Random(), amrex::Random(), amrex::Random()};
    }
};

// struct whose getPositionUnitBox returns x, y and z for a particle with
// regular distribution inside a unit cell.
struct InjectorPositionRegular
{
    InjectorPositionRegular (amrex::Dim3 const& a_ppc) noexcept : ppc(a_ppc) {}

    // i_part: particle number within the cell, required to evenly space
    // particles within the cell.
    // ref_fac: the number of particles evenly-spaced within a cell
    // is a_ppc*(ref_fac**AMREX_SPACEDIM).
    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getPositionUnitBox (int const i_part, int const ref_fac=1) const noexcept
    {
        using namespace amrex;

        int const nx = ref_fac*ppc.x;
        int const ny = ref_fac*ppc.y;
#if (defined WARPX_DIM_3D) || (defined WARPX_DIM_RZ)
        int const nz = ref_fac*ppc.z;
#else
        int const nz = 1;
#endif
        int const ix_part = i_part / (ny*nz);  // written this way backward compatibility
        int const iz_part = (i_part-ix_part*(ny*nz)) / ny;
        int const iy_part = (i_part-ix_part*(ny*nz)) - ny*iz_part;
        return XDim3{
            (0.5_rt + ix_part) / nx,
            (0.5_rt + iy_part) / ny,
            (0.5_rt + iz_part) / nz
        };
    }
private:
    amrex::Dim3 ppc;
};

// Base struct for position injector.
// InjectorPosition contains a union (called Object) that holds any one
// instance of:
// - InjectorPositionRandom : to generate random distribution;
// - InjectorPositionRegular: to generate regular distribution.
// The choice is made at runtime, depending in the constructor called.
// This mimics virtual functions, except the struct is stored in managed memory
// and member functions are made __host__ __device__ to run on CPU and GPU.
// This struct inherits from amrex::Gpu::Managed to provide new and delete
// operators in managed memory when running on GPU. Nothing special on CPU.
struct InjectorPosition
    : public amrex::Gpu::Managed
{
    // This constructor stores a InjectorPositionRandom in union object.
    InjectorPosition (InjectorPositionRandom* t,
                      amrex::Real a_xmin, amrex::Real a_xmax,
                      amrex::Real a_ymin, amrex::Real a_ymax,
                      amrex::Real a_zmin, amrex::Real a_zmax)
        : type(Type::random),
          object(t),
          xmin(a_xmin), xmax(a_xmax),
          ymin(a_ymin), ymax(a_ymax),
          zmin(a_zmin), zmax(a_zmax)
    { }

    // This constructor stores a InjectorPositionRegular in union object.
    InjectorPosition (InjectorPositionRegular* t,
                      amrex::Real a_xmin, amrex::Real a_xmax,
                      amrex::Real a_ymin, amrex::Real a_ymax,
                      amrex::Real a_zmin, amrex::Real a_zmax,
                      amrex::Dim3 const& a_ppc)
        : type(Type::regular),
          object(t, a_ppc),
          xmin(a_xmin), xmax(a_xmax),
          ymin(a_ymin), ymax(a_ymax),
          zmin(a_zmin), zmax(a_zmax)
    { }

    // Explicitly prevent the compiler from generating copy constructors
    // and copy assignment operators.
    InjectorPosition (InjectorPosition const&) = delete;
    InjectorPosition (InjectorPosition&&) = delete;
    void operator= (InjectorPosition const&) = delete;
    void operator= (InjectorPosition &&) = delete;

    // call getPositionUnitBox from the object stored in the union
    // (the union is called Object, and the instance is called object).
    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getPositionUnitBox (int const i_part, int const ref_fac=1) const noexcept
    {
        switch (type)
        {
        case Type::regular:
        {
            return object.regular.getPositionUnitBox(i_part, ref_fac);
        }
        default:
        {
            return object.random.getPositionUnitBox(i_part, ref_fac);
        }
        };
    }

    // bool: whether position specified is within bounds.
    AMREX_GPU_HOST_DEVICE
    bool
    insideBounds (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return (x < xmax and x >= xmin and
                y < ymax and y >= ymin and
                z < zmax and z >= zmin);
    }

private:
    enum struct Type { random, regular };
    Type type;

    // An instance of union Object constructs and stores any one of
    // the objects declared (random or regular).
    union Object {
        Object (InjectorPositionRandom*) noexcept : random() {}
        Object (InjectorPositionRegular*, amrex::Dim3 const& a_ppc) noexcept
            : regular(a_ppc) {}
        InjectorPositionRandom random;
        InjectorPositionRegular regular;
    };
    Object object;

    amrex::Real xmin, xmax;
    amrex::Real ymin, ymax;
    amrex::Real zmin, zmax;
};

#endif
