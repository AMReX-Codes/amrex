#include "QuantumSyncEngineWrapper.H"
//This file provides a wrapper aroud the quantum_sync engine
//provided by the PICSAR library

using namespace picsar::multi_physics;

// Factory class =============================

QuantumSynchrotronEngine::QuantumSynchrotronEngine (){}

QuantumSynchrotronGetOpticalDepth
QuantumSynchrotronEngine::build_optical_depth_functor ()
{
    return QuantumSynchrotronGetOpticalDepth();
}

//============================================
