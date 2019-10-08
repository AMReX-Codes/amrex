#ifndef WARPX_breit_wheeler_engine_wrapper_h_
#define WARPX_breit_wheeler_engine_wrapper_h_

//This file provides a wrapper aroud the breit_wheeler engine
//provided by the QED modules of the PICSAR library

#include "QedWrapperCommons.h"

//BW ENGINE from PICSAR
#include "breit_wheeler_engine.hpp"

using WarpXBreitWheelerWrapper =
    picsar::multi_physics::breit_wheeler_engine<amrex::Real, DummyStruct>;

// Functors ==================================

// These functors provide the core elementary functions of the library
// Can be included in GPU kernels

// Initialization of the optical depth
class BreitWheelerGetOpticalDepth
{
public:
    BreitWheelerGetOpticalDepth()
    {};

    AMREX_GPU_DEVICE
    amrex::Real operator() () const;
};
//____________________________________________

// Factory class =============================
class BreitWheelerEngine
{
public:
    BreitWheelerEngine();

    //Builds the functor to initialize the optical depth
    BreitWheelerGetOpticalDepth build_optical_depth_functor();
};

//============================================

#endif //WARPX_breit_wheeler_engine_wrapper_H_
