#ifndef AMREX_GPU_ASYNC_ARRAY_H_
#define AMREX_GPU_ASYNC_ARRAY_H_

#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <AMReX_Arena.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_GpuDevice.H>

#ifdef AMREX_USE_GPU
extern "C" {
AMREX_HIP_OR_CUDA(
    void  HIPRT_CB amrex_asyncarray_delete ( hipStream_t stream,  hipError_t error, void* p);,
    void CUDART_CB amrex_asyncarray_delete (cudaStream_t stream, cudaError_t error, void* p);)
}
#endif

namespace amrex {
namespace Gpu {

template <typename T, typename = amrex::EnableIf_t<AMREX_IS_TRIVIALLY_COPYABLE(T)> >
class AsyncArray
{
public:

    AsyncArray (T const* h_p, const std::size_t n)
    {
        if (n == 0) return;
        h_data = static_cast<T*>(std::malloc(n*sizeof(T)));
        std::memcpy(h_data, h_p, n*sizeof(T));
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            d_data = static_cast<T*>(The_Device_Arena()->alloc(n*sizeof(T)));
            Gpu::htod_memcpy_async(d_data, h_data, n*sizeof(T));
        }
#endif
    }

    template <typename U = T, typename = amrex::EnableIf_t<std::is_pod<U>::value> >
    explicit AsyncArray (const std::size_t n)
    {
        if (n == 0) return;
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            d_data = static_cast<T*>(The_Device_Arena()->alloc(n*sizeof(T)));
        }
        else
#endif
        {
            h_data = static_cast<T*>(std::malloc(n*sizeof(T)));
        }
    }

    ~AsyncArray () { clear(); }

    T const* data () const noexcept { return (d_data != nullptr) ? d_data : h_data; }
    T* data () noexcept { return (d_data != nullptr) ? d_data : h_data; }
    void clear ()
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            if (d_data != nullptr) {
                T** p = static_cast<T**>(std::malloc(2*sizeof(T*)));
                p[0] = d_data;
                p[1] = h_data;
                AMREX_HIP_OR_CUDA(
                    AMREX_HIP_SAFE_CALL ( hipStreamAddCallback(Gpu::gpuStream(),
                                                               amrex_asyncarray_delete, p, 0));,
                    AMREX_CUDA_SAFE_CALL(cudaStreamAddCallback(Gpu::gpuStream(),
                                                               amrex_asyncarray_delete, p, 0)););
                Gpu::callbackAdded();
            }
        }
        else
#endif
        {
            std::free(h_data);
        }
        d_data = nullptr;
        h_data = nullptr;
    }

    void copyToHost (T* h_p, std::size_t n) const
    {
        if (n == 0) return;
#ifdef AMREX_USE_GPU
        if (d_data)
        {
            Gpu::dtoh_memcpy_async(h_p, d_data, n*sizeof(T));
        }
        else
#endif
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
