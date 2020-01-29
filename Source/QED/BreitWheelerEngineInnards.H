/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
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
    // ctrl contains several parameters:
    // - chi_phot_min : the minium chi parameter to be
    //   considered by the engine
    // - chi_phot_tdndt_min : minimun chi for sub-table 1 (1D)
    // - chi_phot_tdndt_max : maximum chi for sub-table 1 (1D)
    // - chi_phot_tdndt_how_many : how many points to use for sub-table 1 (1D)
    // - chi_phot_tpair_min : minimun chi for sub-table 2 (1D)
    // - chi_phot_tpair_max : maximum chi for sub-table 2 (1D)
    // - chi_phot_tpair_how_many : how many points to use for chi for sub-table 2 (2D)
    // - chi_frac_tpair_how_many : how many points to use for the second axis of sub-table 2 (2D)
    picsar::multi_physics::breit_wheeler_engine_ctrl<amrex::Real> ctrl;

    //Lookup table data
    //---sub-table 1 (1D)
    amrex::Gpu::ManagedVector<amrex::Real> TTfunc_coords;
    amrex::Gpu::ManagedVector<amrex::Real> TTfunc_data;
    //---

    //---sub-table 2 (2D)
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_1;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_2;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_data;
    //______
};
//==========================================================

#endif //WARPX_breit_wheeler_engine_innards_h_
