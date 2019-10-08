#ifndef WARPX_quantum_sync_engine_wrapper_h_
#define WARPX_quantum_sync_engine_wrapper_h_

//This file provides a wrapper aroud the breit_wheeler engine
//provided by the QED modules of the PICSAR library

#include "QedWrapperCommons.h"

//QS ENGINE from PICSAR
#include "quantum_sync_engine.hpp"

using WarpXQuantumSynchrotronWrapper =
    picsar::multi_physics::quantum_synchrotron_engine<amrex::Real, DummyStruct>;

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

// Factory class =============================
class QuantumSynchrotronEngine
{
public:
    QuantumSynchrotronEngine();

    //Builds the functor to initialize the optical depth
    QuantumSynchrotronGetOpticalDepth build_optical_depth_functor();
};

//============================================

#endif //WARPX_quantum_sync_engine_wrapper_h_
