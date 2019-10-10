#include "BreitWheelerEngineWrapper.h"
//This file provides a wrapper aroud the breit_wheeler engine
//provided by the PICSAR library

using namespace picsar::multi_physics;

// Functors ==================================

// Initialization of the optical depth
AMREX_GPU_DEVICE
amrex::Real
BreitWheelerGetOpticalDepth::operator() () const
{
    return WarpXBreitWheelerWrapper::
        internal_get_optical_depth(amrex::Random());
}
//____________________________________________

// Factory class =============================

BreitWheelerEngine::BreitWheelerEngine (){}

BreitWheelerGetOpticalDepth BreitWheelerEngine::build_optical_depth_functor ()
{
    return BreitWheelerGetOpticalDepth();
}

//============================================