#ifndef WARPX_amrex_qed_wrapper_commons_h_
#define WARPX_amrex_qed_wrapper_commons_h_

//Common definitions for the QED library wrappers

#include<AMReX_AmrCore.H>
#include<AMReX_Gpu.H>

#include<utility>

//Sets the decorator for GPU
#define PXRMP_GPU AMREX_GPU_DEVICE
//Sets SI units in the library
#define PXRMP_WITH_SI_UNITS

//An empty data type
struct DummyStruct{};

//An helper function used by the output routines
template <class T>
std::pair<char*, char*> get_begin_end_pointers (T val)
{
    return std::make_pair((char*)&val, ((char*)&val) + sizeof(val));
};

#endif //WARPX_amrex_qed_wrapper_commons_h_
