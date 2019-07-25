#ifndef CUSTOM_MOMENTUM_PROB_H
#define CUSTOM_MOMENTUM_PROB_H

#include <AMReX_ParmParse.H>
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Dim3.H>

// An example of Custom Momentum Profile

struct InjectorMomentumCustom
{
    InjectorMomentumCustom (std::string const& /*a_species_name*/) {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real, amrex::Real, amrex::Real) const noexcept
    {
        return {0., 0., 0.};
    }

    // Note that we are not allowed to have non-trivial destructor.
    // So we rely on clear() to free memory if needed.
    void clear () { }
};

#endif
