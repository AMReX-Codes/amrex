#include "quantum_sync_engine_wrapper.h"
//This file provides a wrapper aroud the quantum_sync engine
//provided by the standard template library

using namespace picsar::multi_physics;

warpx_quantum_sync_engine::warpx_quantum_sync_engine():
    quantum_synchrotron_engine(std::move(amrex_rng_wrapper{}))
{};

//Interface for the get_optical_depth method of the QS engine
amrex::Real
AMREX_GPU_HOST_DEVICE
warpx_quantum_sync_engine::get_optical_depth(){
    return internal_get_optical_depth(amrex::Random());
}
