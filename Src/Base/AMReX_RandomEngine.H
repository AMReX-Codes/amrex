#ifndef AMREX_RANDOM_ENGINE_H_
#define AMREX_RANDOM_ENGINE_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>

#if defined(AMREX_USE_HIP)
#include <hiprand.hpp>
#elif defined(AMREX_USE_CUDA)
#include <curand.h>
#include <curand_kernel.h>
#elif defined(AMREX_USE_DPCPP)
#include <CL/sycl.hpp>
namespace sycl = cl::sycl;
#if (__SYCL_COMPILER_VERSION <= 20200827)
// Disable warnings in mkl's rng header.  Too much noise.
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-braces"
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wextra-semi"
#include <mkl_rng_sycl_device.hpp>
#pragma clang diagnostic pop
#else
#include <oneapi/mkl/rng/device.hpp>
#endif
#if (__SYCL_COMPILER_VERSION >= 20200827)
namespace mkl = oneapi::mkl;
#endif
#endif

namespace amrex
{
#ifdef AMREX_USE_GPU

#ifdef AMREX_USE_DPCPP

    using dpcpp_rng_engine = mkl::rng::device::philox4x32x10<>;
    using dpcpp_rng_descr = mkl::rng::device::engine_descriptor<dpcpp_rng_engine>;
    using dpcpp_rng_acc = mkl::rng::device::engine_accessor<dpcpp_rng_engine>;

    extern dpcpp_rng_descr* rand_engine_descr;

    AMREX_FORCE_INLINE
    dpcpp_rng_descr* getRandEngineDescriptor () { return rand_engine_descr; }

    struct RandomEngine {
        dpcpp_rng_engine* engine;
    };

#else

#ifdef AMREX_USE_HIP
    using randState_t = hiprandState_t;
#else
    using randState_t = curandState_t;
#endif

    extern randState_t* d_states_h_ptr;

    AMREX_FORCE_INLINE
    randState_t* getRandState () { return d_states_h_ptr; }

    struct RandomEngine {
        randState_t* rand_state;
    };

#endif

#else

    struct RandomEngine {}; // CPU

#endif

}

#endif
