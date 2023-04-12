#ifndef AMREX_RANDOM_ENGINE_H_
#define AMREX_RANDOM_ENGINE_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>

#if defined(AMREX_USE_HIP)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#warnings"
#include <hiprand/hiprand.hpp>
#pragma clang diagnostic pop
#elif defined(AMREX_USE_CUDA)
#include <curand.h>
#include <curand_kernel.h>
#elif defined(AMREX_USE_SYCL)
#include <sycl/sycl.hpp>
#include <oneapi/mkl/rng/device.hpp>
namespace mkl = oneapi::mkl;
#endif

namespace amrex
{
#ifdef AMREX_USE_GPU

#ifdef AMREX_USE_SYCL

    using sycl_rng_engine = mkl::rng::device::philox4x32x10<>;
    using sycl_rng_descr = mkl::rng::device::engine_descriptor<sycl_rng_engine>;
    using sycl_rng_acc = mkl::rng::device::engine_accessor<sycl_rng_engine>;

    extern sycl_rng_descr* rand_engine_descr;

    AMREX_FORCE_INLINE
    sycl_rng_descr* getRandEngineDescriptor () { return rand_engine_descr; }

    struct RandomEngine {
        sycl_rng_engine* engine;
    };

#else

#ifdef AMREX_USE_HIP
    using randState_t = hiprandState_t;
#else
    using randState_t = curandState_t;
#endif

    extern randState_t* gpu_rand_state;

    AMREX_FORCE_INLINE
    randState_t* getRandState () { return gpu_rand_state; }

    struct RandomEngine {
        randState_t* rand_state;
    };

#endif

#else

    struct RandomEngine {}; // CPU

#endif

}

#endif
