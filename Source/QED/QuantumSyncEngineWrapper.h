#ifndef WARPX_quantum_sync_engine_wrapper_h_
#define WARPX_quantum_sync_engine_wrapper_h_

#include "QedWrapperCommons.h"

//QS ENGINE from PICSAR
#include "quantum_sync_engine.hpp"

using WarpXQuantumSynchrotronWrapper =
    picsar::multi_physics::quantum_synchrotron_engine<amrex::Real, DummyStruct>;

// Functors ==================================

// These functors provide the core elementary functions of the library
// Can be included in GPU kernels

/* \brief Functor to initialize the optical depth of leptons for the 
*   Quantum Synchrotron process */
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

/* \brief Wrapper for the Quantum Synchrotron engine of the PICSAR library */
class QuantumSynchrotronEngine
{
public:
    QuantumSynchrotronEngine();

    /* \brief Builds the functor to initialize the optical depth */
    QuantumSynchrotronGetOpticalDepth build_optical_depth_functor();
};

//============================================

#endif //WARPX_quantum_sync_engine_wrapper_h_