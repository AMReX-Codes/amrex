#ifndef WARPX_quantum_sync_engine_wrapper_h_
#define WARPX_quantum_sync_engine_wrapper_h_

//This file provides a wrapper aroud the breit_wheeler engine
//provided by the QED modules of the PICSAR library

#include "QedWrapperCommons.h"

//QS ENGINE from PICSAR
#include "quantum_sync_engine.hpp"

using WarpXQuantumSynchrotronWrapper =
    picsar::multi_physics::quantum_synchrotron_engine<amrex::Real, DummyStruct>;

using WarpXQuantumSynchrotronWrapperCtrl =
    picsar::multi_physics::quantum_synchrotron_engine_ctrl<amrex::Real>;

// Struct to hold engine data ================

struct QuantumSynchrotronEngineInnards
{
    // Control parameters
    WarpXQuantumSynchrotronWrapperCtrl ctrl;

    //Lookup table data
    amrex::Gpu::ManagedDeviceVector<amrex::Real> KKfunc_coords;
    amrex::Gpu::ManagedDeviceVector<amrex::Real> KKfunc_data; 
    //______
};

// Functors ==================================

// These functors provide the core elementary functions of the library
// Can be included in GPU kernels

// Initialization of the optical depth
class QuantumSynchrotronGetOpticalDepth
{
public:
    QuantumSynchrotronGetOpticalDepth()
    {};

    AMREX_GPU_DEVICE
    amrex::Real operator() () const;
};
//____________________________________________

// Evolution of the optical depth (returns true if
// an event occurs)
class QuantumSynchrotronEvolveOpticalDepth
{
public:
    QuantumSynchrotronEvolveOpticalDepth(
        QuantumSynchrotronEngineInnards* _innards):
        innards{_innards}{};

    AMREX_GPU_DEVICE
    bool operator()(
    amrex::Real px, amrex::Real py, amrex::Real pz, 
    amrex::Real ex, amrex::Real ey, amrex::Real ez,
    amrex::Real bx, amrex::Real by, amrex::Real bz,
    amrex::Real dt, amrex::Real& opt_depth) const;

private:
    QuantumSynchrotronEngineInnards* innards;
};

// Factory class =============================
class QuantumSynchrotronEngine
{
public:
    QuantumSynchrotronEngine();

    //Builds the functor to initialize the optical depth
    QuantumSynchrotronGetOpticalDepth build_optical_depth_functor();

    //Builds the functor to evolve the optical depth
    QuantumSynchrotronEvolveOpticalDepth build_evolve_functor();

    //Computes the Lookup tables using the default settings 
    //provided by the library
    void computes_lookup_tables_default();

    bool are_lookup_tables_initialized() const;

private:
    bool lookup_tables_initialized = false;

    QuantumSynchrotronEngineInnards innards;

    //Private function which actually computes the lookup tables
    void computes_lookup_tables(
        WarpXQuantumSynchrotronWrapperCtrl ctrl);
};

//============================================

#endif //WARPX_quantum_sync_engine_wrapper_h_
