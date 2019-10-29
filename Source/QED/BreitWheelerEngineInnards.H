#ifndef WARPX_breit_wheeler_engine_innards_h_
#define WARPX_breit_wheeler_engine_innards_h_

#include "QedWrapperCommons.H"

#include <AMReX_Gpu.H>

//This includes only the definition of a simple datastructure
//used to control the Breit Wheeler engine.
#include <breit_wheeler_engine_ctrl.h>

/**
 * This structure holds all the parameters required to use the
 * Breit Wheeler engine: a POD control structure and lookup
 * tables data.
 */
struct BreitWheelerEngineInnards
{
    // Control parameters (a POD struct)
    picsar::multi_physics::breit_wheeler_engine_ctrl<amrex::Real> ctrl;

    //Lookup table data
    amrex::Gpu::ManagedVector<amrex::Real> TTfunc_coords;
    amrex::Gpu::ManagedVector<amrex::Real> TTfunc_data;

    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_1;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_2;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_data;
    //______
};
//==========================================================

#endif //WARPX_breit_wheeler_engine_innards_h_
