#ifndef AMREX_GPU_ASYNC_ARRAY_H_
#define AMREX_GPU_ASYNC_ARRAY_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_GpuAsyncArray.H
 * \brief A GPU array whose device memory is freed asynchronously at stream completion.
 */

#include <AMReX_Arena.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_GpuDevice.H>

#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <memory>

/// \cond DOXYGEN_IGNORE
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
extern "C" {
#if defined(AMREX_USE_HIP)
    void amrex_asyncarray_delete ( hipStream_t stream,  hipError_t error, void* p);
#elif defined(AMREX_USE_CUDA)
    void CUDART_CB amrex_asyncarray_delete (void* p);
#endif
}
#endif
/// \endcond

namespace amrex {
namespace Gpu {

/**
 * \class amrex::Gpu::AsyncArray
 * \ingroup amrex_linear_containers
 * \brief A fixed-size GPU array whose device memory is freed asynchronously
 *        after all pending work on the current GPU stream completes.
 *
 * Useful for short-lived, read-only device arrays.  The host copy lives in
 * pinned memory and is synchronously freed on destruction; the device copy is
 * freed via a stream callback so that any GPU kernels that read it can finish
 * first.
 *
 * \tparam T Element type (must be trivially copyable).
 */
template <typename T, std::enable_if_t<std::is_trivially_copyable_v<T>,int> = 0>
class AsyncArray
{
public:

    //! Construct from a host pointer \p h_p and size \p n; copies data to device asynchronously.
    AsyncArray (T const* h_p, const std::size_t n)
    {
        if (n == 0) { return; }
        h_data = static_cast<T*>(The_Pinned_Arena()->alloc(n*sizeof(T)));
        std::memcpy(h_data, h_p, n*sizeof(T));
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            d_data = static_cast<T*>(The_Arena()->alloc(n*sizeof(T)));
            Gpu::htod_memcpy_async(d_data, h_data, n*sizeof(T));
        }
#endif
    }

    //! Allocate \p n uninitialized elements on the device (or host if not in a launch region).
    template <typename U = T, std::enable_if_t<std::is_standard_layout_v<U> && std::is_trivial_v<U>,int> = 0>
    explicit AsyncArray (const std::size_t n)
    {
        if (n == 0) { return; }
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            d_data = static_cast<T*>(The_Arena()->alloc(n*sizeof(T)));
        }
        else
#endif
        {
            h_data = static_cast<T*>(The_Pinned_Arena()->alloc(n*sizeof(T)));
        }
    }

    //! Destructor: schedules asynchronous deallocation of device memory.
    ~AsyncArray () { clear(); }

    //! AsyncArray is not copyable.
    AsyncArray (AsyncArray const&) = delete;
    //! AsyncArray is not movable.
    AsyncArray (AsyncArray &&) = delete;
    //! AsyncArray is not copy-assignable.
    void operator= (AsyncArray const&) = delete;
    //! AsyncArray is not move-assignable.
    void operator= (AsyncArray &&) = delete;

    //! Return a const pointer to the data (device pointer if available, host otherwise).
    [[nodiscard]] T const* data () const noexcept { return (d_data != nullptr) ? d_data : h_data; }
    //! Return a pointer to the data (device pointer if available, host otherwise).
    [[nodiscard]] T* data () noexcept { return (d_data != nullptr) ? d_data : h_data; }

    //! Asynchronously free device memory and synchronously free host memory.
    void clear ()
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            if (d_data != nullptr) {
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
                T** p = static_cast<T**>(std::malloc(2*sizeof(T*)));
                p[0] = d_data;
                p[1] = h_data;
#if defined(AMREX_USE_HIP)
                AMREX_HIP_SAFE_CALL ( hipStreamAddCallback(Gpu::gpuStream(),
                                                           amrex_asyncarray_delete, p, 0));
#elif defined(AMREX_USE_CUDA)
                AMREX_CUDA_SAFE_CALL(cudaLaunchHostFunc(Gpu::gpuStream(),
                                                        amrex_asyncarray_delete, p));
#endif
#elif defined(AMREX_USE_SYCL)
                auto* pd = d_data;
                auto* ph = h_data;
                auto& q = *(Gpu::gpuStream().queue);
                try {
                    q.submit([&] (sycl::handler& h) {
                        h.host_task([=] () {
                            The_Arena()->free(pd);
                            The_Pinned_Arena()->free(ph);
                        });
                    });
                } catch (sycl::exception const& ex) {
                    amrex::Abort(std::string("host_task: ")+ex.what()+"!!!!!");
                }
#endif
            }
        }
        else
#endif
        {
            The_Pinned_Arena()->free(h_data);
        }
        d_data = nullptr;
        h_data = nullptr;
    }

    //! Copy \p n elements from device memory back to the host pointer \p h_p.
    void copyToHost (T* h_p, std::size_t n) const
    {
        if (n == 0) { return; }
#ifdef AMREX_USE_GPU
        if (d_data)
        {
            Gpu::dtoh_memcpy(h_p, d_data, n*sizeof(T));
        }
        else
#endif
        if (h_data)
        {
            std::memcpy(h_p, h_data, n*sizeof(T));
        }
    }

private:
    T* d_data = nullptr;
    T* h_data = nullptr;
};

}

using Gpu::AsyncArray;
}

#endif
