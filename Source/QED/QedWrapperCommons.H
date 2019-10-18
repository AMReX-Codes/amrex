#ifndef WARPX_amrex_qed_wrapper_commons_h_
#define WARPX_amrex_qed_wrapper_commons_h_

//Common definitions for the QED library wrappers and table builders

#include <AMReX_AmrCore.H>
#include <AMReX_GpuQualifiers.H>

//Sets the decorator for GPU
#define PXRMP_GPU AMREX_GPU_DEVICE
//Sets SI units in the library
#define PXRMP_WITH_SI_UNITS

namespace QedUtils{
    //An empty data type
    struct DummyStruct{};
};

#endif //WARPX_amrex_qed_wrapper_commons_h_
