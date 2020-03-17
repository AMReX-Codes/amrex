/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_quantum_sync_engine_innards_h_
#define WARPX_quantum_sync_engine_innards_h_

#include "QedWrapperCommons.H"

#include <AMReX_Gpu.H>

//This includes only the definition of a simple datastructure
//used to control the Quantum Synchrotron engine.
#include <quantum_sync_engine_ctrl.h>

/**
 * This structure holds all the parameters required to use the
 * Quantum Synchrotron engine: a POD control structure and lookup
 * tables data.
 */
struct QuantumSynchrotronEngineInnards
{
    // Control parameters (a POD struct)
    // ctrl contains several parameters:
    // - chi_part_min : the minium chi parameter to be
    //   considered by the engine
    // - chi_part_tdndt_min : minimun chi for sub-table 1 (1D)
    // - chi_part_tdndt_max : maximum chi for sub-table 1 (1D)
    // - chi_part_tdndt_how_many : how many points to use for sub-table 1 (1D)
    // - chi_part_tem_min : minimun chi for sub-table 2 (1D)
    // - chi_part_tem_max : maximum chi for sub-table 2 (1D)
    // - chi_part_tem_how_many : how many points to use for chi for sub-table 2 (2D)
    // - prob_tem_how_many : how many points to use for the second axis of sub-table 2 (2D)
    picsar::multi_physics::quantum_synchrotron_engine_ctrl<amrex::Real> ctrl;

    //Lookup table data
    //---sub-table 1 (1D)
    amrex::Gpu::ManagedDeviceVector<amrex::Real> KKfunc_coords;
    amrex::Gpu::ManagedDeviceVector<amrex::Real> KKfunc_data;
    //---

    //---sub-table 2 (2D)
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_1;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_2;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_data;
    //______
};
//==========================================================

#endif //WARPX_quantum_sync_engine_innards_h_
