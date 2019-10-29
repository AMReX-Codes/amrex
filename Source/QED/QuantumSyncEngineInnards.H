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
    picsar::multi_physics::quantum_synchrotron_engine_ctrl<amrex::Real> ctrl;

    //Lookup table data (should stay on GPU)
    amrex::Gpu::ManagedVector<amrex::Real> KKfunc_coords;
    amrex::Gpu::ManagedVector<amrex::Real> KKfunc_data;

    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_1;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_2;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_data;
    //______
};
//==========================================================

#endif //WARPX_quantum_sync_engine_innards_h_
