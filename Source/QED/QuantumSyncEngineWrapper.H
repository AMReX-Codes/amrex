#ifndef WARPX_quantum_sync_engine_wrapper_h_
#define WARPX_quantum_sync_engine_wrapper_h_

#include "QedWrapperCommons.H"

//QS ENGINE from PICSAR
#include "quantum_sync_engine.hpp"

using WarpXQuantumSynchrotronWrapper =
    picsar::multi_physics::quantum_synchrotron_engine<amrex::Real, DummyStruct>;

// Functors ==================================

// These functors provide the core elementary functions of the library
// Can be included in GPU kernels

/**
*  Functor to initialize the optical depth of leptons for the
*  Quantum Synchrotron process
*/
class QuantumSynchrotronGetOpticalDepth
{
public:
    QuantumSynchrotronGetOpticalDepth ()
    {};

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real operator() () const
    {
        return WarpXQuantumSynchrotronWrapper::
            internal_get_optical_depth(amrex::Random());
    }
};
//____________________________________________

// Factory class =============================

/**
 * \brief Wrapper for the Quantum Synchrotron engine of the PICSAR library
 */
class QuantumSynchrotronEngine
{
public:
    QuantumSynchrotronEngine ();

    /**
     * \brief Builds the functor to initialize the optical depth
     */
    QuantumSynchrotronGetOpticalDepth build_optical_depth_functor ();
};

//============================================

#endif //WARPX_quantum_sync_engine_wrapper_h_
