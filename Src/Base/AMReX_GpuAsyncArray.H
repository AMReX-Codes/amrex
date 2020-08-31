#ifndef AMREX_GPU_ASYNC_ARRAY_H_
#define AMREX_GPU_ASYNC_ARRAY_H_

#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <AMReX_Arena.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_GpuDevice.H>

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
extern "C" {
#if defined(AMREX_USE_HIP)
    void amrex_asyncarray_delete ( hipStream_t stream,  hipError_t error, void* p);
#elif defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10)
    void CUDART_CB amrex_asyncarray_delete (void* p);
#elif defined(AMREX_USE_CUDA)
    void CUDART_CB amrex_asyncarray_delete (cudaStream_t stream, cudaError_t error, void* p);
#endif
}
#endif

namespace amrex {
namespace Gpu {

template <typename T, typename std::enable_if<AMREX_IS_TRIVIALLY_COPYABLE(T),int>::type = 0>
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

    template <typename U = T, typename std::enable_if<std::is_pod<U>::value,int>::type = 0>
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
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
                T** p = static_cast<T**>(std::malloc(2*sizeof(T*)));
                p[0] = d_data;
                p[1] = h_data;
#if defined(AMREX_USE_HIP)
                AMREX_HIP_SAFE_CALL ( hipStreamAddCallback(Gpu::gpuStream(),
                                                           amrex_asyncarray_delete, p, 0));
#elif defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10)
                AMREX_CUDA_SAFE_CALL(cudaLaunchHostFunc(Gpu::gpuStream(),
                                                        amrex_asyncarray_delete, p));
#elif defined(AMREX_USE_CUDA)
                AMREX_CUDA_SAFE_CALL(cudaStreamAddCallback(Gpu::gpuStream(),
                                                           amrex_asyncarray_delete, p, 0));
#endif
                Gpu::callbackAdded();
#else
                // xxxxx DPCPP todo
                Gpu::streamSynchronize();
                The_Device_Arena()->free(d_data);
                std::free(h_data);
#endif
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
