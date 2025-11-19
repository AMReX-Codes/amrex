#ifndef AMREX_RANDOM_ENGINE_H_
#define AMREX_RANDOM_ENGINE_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#ifdef AMREX_USE_OMP
#  include <AMReX_BLassert.H>
#  include <omp.h>
#endif

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
#include <oneapi/mkl/rng.hpp>
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

    inline
    sycl_rng_descr* getRandEngineDescriptor () {
#ifdef AMREX_USE_OMP
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(! omp_in_parallel(), "RNG on GPU is not OpenMP thread safe");
#endif
        return rand_engine_descr;
    }

    struct RandomEngine {
        sycl_rng_engine* engine;
        AMREX_GPU_HOST_DEVICE RandomEngine (sycl_rng_engine* e) : engine(e) {}
    };

#else

#ifdef AMREX_USE_HIP
    using randState_t = hiprandState_t;
    using randGenerator_t = hiprandGenerator_t;
#else
    using randState_t = curandState_t;
    using randGenerator_t = curandGenerator_t;
#endif

    extern randState_t* gpu_rand_state;

    inline
    randState_t* getRandState () {
#ifdef AMREX_USE_OMP
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(! omp_in_parallel(), "RNG on GPU is not OpenMP thread safe");
#endif
        return gpu_rand_state;
    }

    struct RandomEngine {
        randState_t* rand_state;
        AMREX_GPU_HOST_DEVICE RandomEngine (randState_t* rs) : rand_state(rs) {}
    };

#endif

#else

    struct RandomEngine {}; // CPU

#endif

    AMREX_FORCE_INLINE
    RandomEngine getInvalidRandomEngine () {
#ifdef AMREX_USE_GPU
        return RandomEngine{nullptr};
#else
        return RandomEngine{};
#endif
    }

}

#endif
