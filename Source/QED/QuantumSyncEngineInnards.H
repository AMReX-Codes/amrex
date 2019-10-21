#ifndef WARPX_quantum_sync_engine_innards_h_
#define WARPX_quantum_sync_engine_innards_h_

#include "QedWrapperCommons.H"

#include <AMReX_Gpu.H>

#include <quantum_sync_engine_ctrl.h>

// Struct to hold Quantum Synchrotron engine data ================

struct QuantumSynchrotronEngineInnards
{
    // Control parameters
    picsar::multi_physics::quantum_synchrotron_engine_ctrl<amrex::Real> ctrl;

    //Lookup table data
    amrex::Gpu::ManagedDeviceVector<amrex::Real> KKfunc_coords;
    amrex::Gpu::ManagedDeviceVector<amrex::Real> KKfunc_data;

    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_1;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_2;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_data;
    //______
};


#endif //WARPX_quantum_sync_engine_innards_h_