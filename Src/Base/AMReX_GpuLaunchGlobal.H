#ifndef AMREX_GPU_LAUNCH_GLOBAL_H_
#define AMREX_GPU_LAUNCH_GLOBAL_H_

#include <AMReX_GpuQualifiers.H>

namespace amrex {

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    // Because of CUDA, we cannot take rvalue lambdas.
    template<class L>
    AMREX_GPU_GLOBAL void launch_global (L f0) { f0(); }
#endif

}

#endif
