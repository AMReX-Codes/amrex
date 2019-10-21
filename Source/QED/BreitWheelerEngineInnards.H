#ifndef WARPX_breit_wheeler_engine_innards_h_
#define WARPX_breit_wheeler_engine_innards_h_

#include "QedWrapperCommons.H"

#include <AMReX_Gpu.H>

#include <breit_wheeler_engine_ctrl.h>

// Struct to hold Breit Wheeler engine data ================
struct BreitWheelerEngineInnards
{
    // Control parameters
    picsar::multi_physics::breit_wheeler_engine_ctrl<amrex::Real> ctrl;

    //Lookup table data
    amrex::Gpu::ManagedVector<amrex::Real> TTfunc_coords;
    amrex::Gpu::ManagedVector<amrex::Real> TTfunc_data;

    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_1;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_2;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_data;
    //______
};

#endif //WARPX_breit_wheeler_engine_innards_h_