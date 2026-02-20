#ifndef AMREX_GPU_TYPES_H_
#define AMREX_GPU_TYPES_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_GpuTypes.H
 * \brief GPU execution types: dim3, gpuStream_t, and the Gpu::Handler struct.
 */

#include <AMReX_Extension.H>

#ifdef AMREX_USE_GPU

#ifdef AMREX_USE_SYCL
#  include <sycl/sycl.hpp>
#endif

namespace amrex {

#ifdef AMREX_USE_SYCL

/// \ingroup amrex_execution
/// 3D unsigned integer grid/block dimension (SYCL equivalent of CUDA dim3).
struct dim3 {
    unsigned int x = 1; ///< x dimension.
    unsigned int y = 1; ///< y dimension.
    unsigned int z = 1; ///< z dimension.
    //! Default-construct with (1,1,1).
    dim3 () = default;
    //! Construct with given x, y, z dimensions.
    constexpr dim3 (unsigned int x_, unsigned int y_=1, unsigned int z_=1) : x(x_),y(y_),z(z_) {}
};

/// \ingroup amrex_execution
/// 1D grid dimension used internally for SYCL kernel launches.
struct Dim1 {
    std::size_t x; ///< x dimension.
};

/// \ingroup amrex_execution
/// SYCL equivalent of a CUDA/HIP stream handle.
struct gpuStream_t {
    sycl::queue* queue = nullptr; ///< SYCL queue handle.
    //! Return true if both handles refer to the same SYCL queue.
    bool operator== (gpuStream_t const& rhs) const noexcept { return queue == rhs.queue; }
    //! Return true if both handles refer to different SYCL queues.
    bool operator!= (gpuStream_t const& rhs) const noexcept { return queue != rhs.queue; }
};

#endif

}

#endif

namespace amrex::Gpu {

/**
 * \ingroup amrex_execution
 * \brief Kernel handler passed to GPU kernels to expose thread/block indices
 *        and synchronization primitives in a backend-portable way.
 *
 * In CUDA/HIP builds the hardware intrinsics (threadIdx, blockIdx, etc.) are
 * used directly inside a kernel, so this struct only carries the number of
 * active threads when the block is not full. In SYCL builds it wraps the
 * \c sycl::nd_item and shared local memory pointer.
 */
#if defined(AMREX_USE_SYCL)

struct Handler
{
    AMREX_GPU_HOST_DEVICE constexpr
    Handler (sycl::nd_item<1> const* a_item = nullptr, void* a_local = nullptr,
             int a_n_active_threds = -1)
        : item(a_item), local(a_local), numActiveThreads(a_n_active_threds) {}

    //! Return true if all threads in the work-group are active.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isFullBlock () const {
        return (numActiveThreads >= int(item->get_local_range(0)) ||
                numActiveThreads <= 0);
    }

    //! Return the global (flattened) thread index.
    std::size_t globalIdx () const { return item->get_global_linear_id(); }
    //! Return the work-group (block) index.
    std::size_t blockIdx () const { return item->get_group_linear_id(); }
    //! Return the local (thread-within-block) index.
    std::size_t threadIdx () const { return item->get_local_linear_id(); }
    //! Return the number of work-groups (grid dimension).
    std::size_t gridDim () const { return item->get_group_range(0); }
    //! Return the work-group size (block dimension).
    std::size_t blockDim () const { return item->get_local_range(0); }
    //! Return the sub-group (warp) index within the work-group.
    std::size_t warpIdx () const { return item->get_sub_group().get_group_id()[0]; }
    //! Return the lane index within the sub-group (warp).
    std::size_t laneIdx () const { return item->get_sub_group().get_local_id()[0]; }
    //! Return the sub-group (warp) size.
    std::size_t warpDim () const { return item->get_sub_group().get_group_range()[0]; }

    //! Return a pointer to the shared local memory.
    void* sharedMemory () const { return local; }

    //! Synchronize across shared local memory.
    void sharedBarrier () const { item->barrier(sycl::access::fence_space::local_space); }
    //! Synchronize across global memory.
    void globalBarrier () const { item->barrier(sycl::access::fence_space::global_space); }
    //! Synchronize across both local and global memory.
    void syncThreads () const { item->barrier(sycl::access::fence_space::global_and_local); }

    sycl::nd_item<1> const* item;
    void* local; ///< SYCL shared local memory pointer.
    int numActiveThreads; ///< Number of active threads; <= 0 means all threads are active.
};

#elif defined(AMREX_USE_GPU)

struct Handler
{
    //! Construct with an optional active-thread count (\p n_active_threads <= 0 means all active).
    AMREX_GPU_HOST_DEVICE constexpr explicit Handler (int n_active_threads = -1)
        : numActiveThreads(n_active_threads) {}

    //! Return true if all threads in the CUDA/HIP block are active.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isFullBlock () const {
        return (numActiveThreads >= (int)blockDim.x ||
                numActiveThreads <= 0);
    }

    int numActiveThreads; ///< Number of active threads; <= 0 means all threads are active.
};

#else

/// \ingroup amrex_execution
//! CPU stub for Handler: no GPU state to expose.
struct Handler {};

#endif

}

#endif
