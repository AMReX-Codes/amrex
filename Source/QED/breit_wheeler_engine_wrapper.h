#ifndef WARPX_breit_wheeler_engine_wrapper_h_
#define WARPX_breit_wheeler_engine_wrapper_h_

//This file provides a wrapper aroud the breit_wheeler engine
//provided by the standard template library

//BW ENGINE
//#define PXRMP_GPU __host__ __device__
#define PXRMP_WITH_SI_UNITS
#include "breit_wheeler_engine.hpp"

#include "amrex_rng_wrapper.h"

using warpx_breit_wheeler_engine =
  picsar::multi_physics::breit_wheeler_engine<amrex::Real, amrex_rng_wrapper>;

//Helper function to initialize the engine
inline warpx_breit_wheeler_engine init_warpx_breit_wheeler_engine(){
  return  warpx_breit_wheeler_engine{std::move(amrex_rng_wrapper{})};
}

//Interface for the get_optical_depth method of the BW engine
inline
AMREX_GPU_HOST_DEVICE
amrex::Real warpx_breit_wheeler_get_optical_depth(){
    return warpx_breit_wheeler_engine::
        internal_get_optical_depth(amrex::Random());
}


//___________________________________________

#endif //WARPX_breit_wheeler_engine_wrapper_H_
