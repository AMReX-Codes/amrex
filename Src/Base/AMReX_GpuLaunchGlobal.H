#ifndef AMREX_GPU_LAUNCH_GLOBAL_H_
#define AMREX_GPU_LAUNCH_GLOBAL_H_

#include <AMReX_GpuQualifiers.H>

namespace amrex {

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    template<class L>
    AMREX_GPU_GLOBAL void launch_global (L f0) { f0(); }

    template<class L, class... Lambdas>
    AMREX_GPU_GLOBAL void launch_global (L f0, Lambdas... fs) { f0(); call_device(fs...); }
#endif

}

#endif
