/* Copyright 2019 Maxence Thevenet, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef CUSTOM_MOMENTUM_PROB_H
#define CUSTOM_MOMENTUM_PROB_H

#include <AMReX_ParmParse.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Dim3.H>

// An example of Custom Momentum Profile

// struct whose getDensity returns momentum at a given position computed from
// a custom function.
struct InjectorMomentumCustom
{
    InjectorMomentumCustom (std::string const& /*a_species_name*/) {}

    // Return momentum at given position (illustration: momentum=0).
    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real, amrex::Real, amrex::Real) const noexcept
    {
        return {0., 0., 0.};
    }

    // Return momentum at given position (illustration: momentum=0).
    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getBulkMomentum (amrex::Real, amrex::Real, amrex::Real) const noexcept
    {
        return {0., 0., 0.};
    }

    // Note that we are not allowed to have non-trivial destructor.
    // So we rely on clear() to free memory if needed.
    void clear () { }
};

#endif
