#ifndef WARPX_breit_wheeler_engine_wrapper_h_
#define WARPX_breit_wheeler_engine_wrapper_h_

#include "QedWrapperCommons.H"

//BW ENGINE from PICSAR
#include "breit_wheeler_engine.hpp"

using WarpXBreitWheelerWrapper =
    picsar::multi_physics::breit_wheeler_engine<amrex::Real, DummyStruct>;

// Functors ==================================

// These functors provide the core elementary functions of the library
// Can be included in GPU kernels

/**
 * \brief Functor to initialize the optical depth of photons for the
 * Breit-Wheeler process
 */
class BreitWheelerGetOpticalDepth
{
public:
    BreitWheelerGetOpticalDepth ()
    {};

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real operator() () const
    {
        return WarpXBreitWheelerWrapper::
            internal_get_optical_depth(amrex::Random());
    }
};
//____________________________________________

// Factory class =============================

/**
 * \brief Wrapper for the Breit Wheeler engine of the PICSAR library
 */
class BreitWheelerEngine
{
public:
    BreitWheelerEngine ();

    /**
     * \brief Builds the functor to initialize the optical depth
     */
    BreitWheelerGetOpticalDepth build_optical_depth_functor ();
};

//============================================

#endif //WARPX_breit_wheeler_engine_wrapper_H_
