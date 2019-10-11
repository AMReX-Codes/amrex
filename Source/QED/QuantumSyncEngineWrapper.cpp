#include "QuantumSyncEngineWrapper.hpp"
//This file provides a wrapper aroud the quantum_sync engine
//provided by the PICSAR library

using namespace picsar::multi_physics;

// Functors ==================================

// Initialization of the optical depth

AMREX_GPU_DEVICE
amrex::Real
QuantumSynchrotronGetOpticalDepth::operator() () const
{
    return WarpXQuantumSynchrotronWrapper::
        internal_get_optical_depth(amrex::Random());
}
//____________________________________________

// Factory class =============================

QuantumSynchrotronEngine::QuantumSynchrotronEngine (){}

QuantumSynchrotronGetOpticalDepth
QuantumSynchrotronEngine::build_optical_depth_functor ()
{
    return QuantumSynchrotronGetOpticalDepth();
}

//============================================