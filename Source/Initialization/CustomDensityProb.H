/* Copyright 2019 Maxence Thevenet, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef CUSTOM_DENSITY_PROB_H_
#define CUSTOM_DENSITY_PROB_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>
#include <AMReX_Dim3.H>

// An example of Custom Density Profile

// struct whose getDensity returns density at a given position computed from
// a custom function, with runtime input parameters.
struct InjectorDensityCustom
{
    InjectorDensityCustom (std::string const& species_name)
        : p(nullptr)
    {
        // Read parameters for custom density profile from file, and
        // store them in managed memory.
        amrex::ParmParse pp(species_name);
        std::vector<amrex::Real> v;
        pp.getarr("custom_profile_params", v);
        p = static_cast<amrex::Real*>
            (amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real)*v.size()));
        for (int i = 0; i < static_cast<int>(v.size()); ++i) {
            p[i] = v[i];
        }
    }

    // Return density at given position, using user-defined parameters
    // stored in p.
    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getDensity (amrex::Real, amrex::Real, amrex::Real) const noexcept
    {
        return p[0];
    }

    // Note that we are not allowed to have non-trivial destructor.
    // So we rely on clear() to free memory.
    void clear () {
        amrex::The_Managed_Arena()->free(p);
    }

private:
    amrex::Real* p;
};

#endif
