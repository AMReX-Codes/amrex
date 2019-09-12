#ifndef WARPX_quantum_sync_engine_wrapper_h_
#define WARPX_quantum_sync_engine_wrapper_h_

//This file provides a wrapper aroud the quantum_synchrotron_engine
//provided by the standard template library

//BW ENGINE
#include "qed_wrapper_commons.h"
#include "quantum_sync_engine.hpp"

#include "amrex_rng_wrapper.h"

class warpx_quantum_sync_engine :
    public picsar::multi_physics::
    quantum_synchrotron_engine<amrex::Real, amrex_rng_wrapper>
{
public:
    warpx_quantum_sync_engine();

    //Interface for the get_optical_depth method of the BW engine
    static AMREX_GPU_HOST_DEVICE
    amrex::Real get_optical_depth();
};

//___________________________________________

#endif //WARPX_quantum_sync_engine_wrapper_h_
