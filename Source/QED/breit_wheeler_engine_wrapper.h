#ifndef WARPX_breit_wheeler_engine_wrapper_h_
#define WARPX_breit_wheeler_engine_wrapper_h_

//This file provides a wrapper aroud the breit_wheeler engine
//provided by the standard template library

//BW ENGINE
//#define PXRMP_GPU __host__ __device__
#define PXRMP_WITH_SI_UNITS
#include "breit_wheeler_engine.hpp"

#include "amrex_rng_wrapper.hpp"

using warpx_breit_wheeler_engine =
  picsar::multi_physics::breit_wheeler_engine<amrex::Real, amrex_rng_wrapper>;

//Helper function to initialize the engine
inline warpx_breit_wheeler_engine&& init_warpx_breit_wheeler_engine(){
  return std::move(warpx_breit_wheeler_engine{std::move{amrex_rng_wrapper{}}});
}

//___________________________________________
#endif

#endif //WARPX_breit_wheeler_engine_wrapper_H_
