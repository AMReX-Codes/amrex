#include "BreitWheelerEngineWrapper.H"
//This file provides a wrapper aroud the breit_wheeler engine
//provided by the PICSAR library

using namespace picsar::multi_physics;

// Factory class =============================

BreitWheelerEngine::BreitWheelerEngine (){}

BreitWheelerGetOpticalDepth BreitWheelerEngine::build_optical_depth_functor ()
{
    return BreitWheelerGetOpticalDepth();
}

//============================================
