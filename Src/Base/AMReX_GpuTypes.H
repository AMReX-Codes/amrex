#ifndef AMREX_GPU_TYPES_H_
#define AMREX_GPU_TYPES_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>

#ifdef AMREX_USE_GPU

#ifdef AMREX_USE_SYCL
#  include <sycl/sycl.hpp>
#endif

namespace amrex {

#ifdef AMREX_USE_SYCL

struct dim3 {
    unsigned int x = 1;
    unsigned int y = 1;
    unsigned int z = 1;
    dim3 () = default;
    constexpr dim3 (unsigned int x_, unsigned int y_=1, unsigned int z_=1) : x(x_),y(y_),z(z_) {}
};

struct Dim1 {
    std::size_t x;
};

struct gpuStream_t {
    sycl::queue* queue = nullptr;
    bool operator== (gpuStream_t const& rhs) noexcept { return queue == rhs.queue; }
};

#endif

}

#endif

namespace amrex::Gpu {

#if defined(AMREX_USE_SYCL)

struct Handler
{
    AMREX_GPU_HOST_DEVICE constexpr
    Handler (sycl::nd_item<1> const* a_item = nullptr, void* a_local = nullptr,
             int a_n_active_threds = -1)
        : item(a_item), local(a_local), numActiveThreads(a_n_active_threds) {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isFullBlock () const {
        return (numActiveThreads >= int(item->get_local_range(0)) ||
                numActiveThreads <= 0);
    }

    std::size_t globalIdx () const { return item->get_global_linear_id(); }
    std::size_t blockIdx () const { return item->get_group_linear_id(); }
    std::size_t threadIdx () const { return item->get_local_linear_id(); }
    //
    std::size_t gridDim () const { return item->get_group_range(0); }
    std::size_t blockDim () const { return item->get_local_range(0); }

    // warp index in block
    std::size_t warpIdx () const { return item->get_sub_group().get_group_id()[0]; }
    // lane index in warp
    std::size_t laneIdx () const { return item->get_sub_group().get_local_id()[0]; }
    // warp size
    std::size_t warpDim () const { return item->get_sub_group().get_group_range()[0]; }

    void* sharedMemory () const { return local; }

    void sharedBarrier () const { item->barrier(sycl::access::fence_space::local_space); }
    void globalBarrier () const { item->barrier(sycl::access::fence_space::global_space); }
    void syncThreads () const { item->barrier(sycl::access::fence_space::global_and_local); }

    sycl::nd_item<1> const* item;
    void* local; // SYCL shared local memory
    int numActiveThreads;
};

#elif defined(AMREX_USE_GPU)

struct Handler
{
    AMREX_GPU_HOST_DEVICE constexpr Handler (int n_active_threads = -1)
        : numActiveThreads(n_active_threads) {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isFullBlock () const {
        return (numActiveThreads >= (int)blockDim.x ||
                numActiveThreads <= 0);
    }

    int numActiveThreads;
};

#else

struct Handler {};

#endif

}

#endif
