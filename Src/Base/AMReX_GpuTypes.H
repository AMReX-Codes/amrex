#ifndef AMREX_GPU_TYPES_H_
#define AMREX_GPU_TYPES_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_BaseFwd.H>
#include <AMReX_GpuControl.H>

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
    bool operator!= (gpuStream_t const& rhs) noexcept { return queue != rhs.queue; }
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
    AMREX_GPU_HOST_DEVICE constexpr explicit Handler (int n_active_threads = -1)
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

template<unsigned int threads_per_block, typename T, int dim>
struct LaunchHandler
{
    static_assert(dim >= 1 && dim <= 3);

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    LaunchHandler() = default;
#elif defined(AMREX_USE_SYCL)
    LaunchHandler(sycl::nd_item<dim> const* a_item, T * a_shared_mem)
        : m_item{a_item}, m_shared_mem{a_shared_mem} {}
#else
    LaunchHandler(IntVectND<dim> a_blockid, IntVectND<dim> a_griddim, T * a_shared_mem)
        : m_blockid{a_blockid}, m_griddim{a_griddim}, m_shared_mem{a_shared_mem} {}
#endif

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    unsigned int threadIdx1D () const {
        return static_cast<unsigned int>(AMREX_HIP_CUDA_OR_SYCL_OR_CPU(
            threadIdx.x,
            m_item->get_local_linear_id(),
            0
        ));
    }

    template<int ...nd_block_size>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> threadIdxND () const {
        static_assert(sizeof...(nd_block_size) == dim &&
                      (1 * ... * nd_block_size) == threads_per_block);
        constexpr IntVectND<dim> iv_block_size = blockDimND<nd_block_size...>();
        IntVectND<dim> ret(0);
        unsigned int idx = threadIdx1D();
        if constexpr (dim == 3) {
            ret[2] = idx / iv_block_size[1];
            idx = idx - ret[2] * iv_block_size[1];
        }
        if constexpr (dim >= 2) {
            ret[1] = idx / iv_block_size[0];
            idx = idx - ret[1] * iv_block_size[0];
        }
        ret[0] = idx;
        return ret;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    unsigned int blockIdx1D () const {
        if constexpr (dim == 1) {
            return static_cast<unsigned int>(AMREX_HIP_CUDA_OR_SYCL_OR_CPU(
                blockIdx.x,
                m_item->get_group_id(0),
                m_blockid[0]
            ));
        } else if constexpr (dim == 2) {
            return static_cast<unsigned int>(AMREX_HIP_CUDA_OR_SYCL_OR_CPU(
                blockIdx.x + gridDim.x * blockIdx.y,
                m_item->get_group_linear_id(),
                m_blockid[0] + m_griddim[0] * m_blockid[1]
            ));
        } else {
            return static_cast<unsigned int>(AMREX_HIP_CUDA_OR_SYCL_OR_CPU(
                blockIdx.x + gridDim.x * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z,
                m_item->get_group_linear_id(),
                m_blockid[0] + m_griddim[0] * m_blockid[1] +
                m_griddim[0] * m_griddim[1] * m_blockid[2]
            ));
        }
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> blockIdxND () const {
        if constexpr (dim == 1) {
            return AMREX_HIP_CUDA_OR_SYCL_OR_CPU(
                IntVectND<dim>(static_cast<int>(blockIdx.x)),
                IntVectND<dim>(static_cast<int>(m_item->get_group_id(0))),
                m_blockid
            );
        } else if constexpr (dim == 2) {
            return AMREX_HIP_CUDA_OR_SYCL_OR_CPU(
                IntVectND<dim>(static_cast<int>(blockIdx.x),
                               static_cast<int>(blockIdx.y)),
                IntVectND<dim>(static_cast<int>(m_item->get_group_id(1)),
                               static_cast<int>(m_item->get_group_id(0))),
                m_blockid
            );
        } else {
            return AMREX_HIP_CUDA_OR_SYCL_OR_CPU(
                IntVectND<dim>(static_cast<int>(blockIdx.x),
                               static_cast<int>(blockIdx.y),
                               static_cast<int>(blockIdx.z)),
                IntVectND<dim>(static_cast<int>(m_item->get_group_id(2)),
                               static_cast<int>(m_item->get_group_id(1)),
                               static_cast<int>(m_item->get_group_id(0))),
                m_blockid
            );
        }
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static constexpr unsigned int blockDim1D () {
        return threads_per_block;
    }

    template<int ...nd_block_size>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static constexpr IntVectND<dim> blockDimND () {
        static_assert(sizeof...(nd_block_size) == dim &&
                      (1 * ... * nd_block_size) == threads_per_block);
        return {nd_block_size...};
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    unsigned int globalIdx1D () const {
        return blockIdx1D() * threads_per_block + threadIdx1D();
    }

    template<int ...nd_block_size>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> globalIdxND () const {
        return blockIdxND() * blockDimND<nd_block_size...>() + threadIdxND<nd_block_size...>();
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void syncthreads () const {
        AMREX_HIP_CUDA_OR_SYCL_OR_CPU(
            __syncthreads(),
            m_item->barrier(sycl::access::fence_space::global_and_local),
            (void)0
        );
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T* shared_memory () const {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        static_assert(32 >= alignof(T));
        alignas(32) extern __shared__ char smem[];
        return reinterpret_cast<T*>(smem);
#else
        return m_shared_mem;
#endif
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto handler () const {
#if defined(AMREX_USE_SYCL)
        return Gpu::Handler(m_item, reinterpret_cast<void*>(m_shared_mem), threads_per_block);
#elif defined(AMREX_USE_GPU)
        return Gpu::Handler(threads_per_block);
#else
        return 0;
#endif
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto sycl_item () const {
#if defined(AMREX_USE_SYCL)
        return m_item;
#else
        return 0;
#endif
    }

private:

#if defined(AMREX_USE_SYCL)
    sycl::nd_item<dim> const* m_item;
    T * m_shared_mem;
#elif !defined(AMREX_USE_GPU)
    IntVectND<dim> m_blockid;
    IntVectND<dim> m_griddim;
    T * m_shared_mem;
#endif
};

}

#endif
