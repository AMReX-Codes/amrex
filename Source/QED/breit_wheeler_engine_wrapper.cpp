#include "breit_wheeler_engine_wrapper.h"
//This file provides a wrapper aroud the breit_wheeler engine
//provided by the standard template library

using namespace picsar::multi_physics;

warpx_breit_wheeler_engine::warpx_breit_wheeler_engine():
    breit_wheeler_engine(std::move(amrex_rng_wrapper{}))
{};

//Interface for the get_optical_depth method of the BW engine
amrex::Real
AMREX_GPU_HOST_DEVICE
warpx_breit_wheeler_engine::get_optical_depth(){
    return internal_get_optical_depth(amrex::Random());
}
