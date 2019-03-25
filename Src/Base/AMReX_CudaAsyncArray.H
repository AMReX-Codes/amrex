#ifndef AMREX_CUDA_ASYNC_ARRAY_H_
#define AMREX_CUDA_ASYNC_ARRAY_H_

#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <AMReX_Arena.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_CudaDevice.H>

#ifdef AMREX_USE_CUDA
extern "C" {
    void CUDART_CB amrex_asyncarray_delete (cudaStream_t stream, cudaError_t error, void* p);
}
#endif

namespace amrex {
namespace Cuda {

template <typename T, typename = amrex::EnableIf_t<AMREX_IS_TRIVIALLY_COPYABLE(T)> >
class AsyncArray
{
public:

    AsyncArray (T const* h_p, const std::size_t n)
    {
        if (n == 0) return;
        h_data = static_cast<T*>(std::malloc(n*sizeof(T)));
        std::memcpy(h_data, h_p, n*sizeof(T));
#ifdef AMREX_USE_CUDA
        if (inLaunchRegion())
        {
            d_data = static_cast<T*>(The_Device_Arena()->alloc(n*sizeof(T)));
            AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(d_data, h_data, n*sizeof(T),
                                                cudaMemcpyHostToDevice, Device::cudaStream()));
        }
#endif
    }

    template <typename U = T, typename = amrex::EnableIf_t<std::is_pod<U>::value> >
    AsyncArray (const std::size_t n)
    {
        if (n == 0) return;
#ifdef AMREX_USE_CUDA
        if (inLaunchRegion())
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
#ifdef AMREX_USE_CUDA
        if (inLaunchRegion())
        {
            if (d_data != nullptr) {
                T** p = static_cast<T**>(std::malloc(2*sizeof(T*)));
                p[0] = d_data;
                p[1] = h_data;
                AMREX_GPU_SAFE_CALL(cudaStreamAddCallback(Device::cudaStream(),
                                                          amrex_asyncarray_delete, p, 0));
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
#ifdef AMREX_USE_CUDA
        if (d_data)
        {
            AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(h_p, d_data, n*sizeof(T),
                                                cudaMemcpyDeviceToHost, Device::cudaStream()));
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

using Cuda::AsyncArray;
}

#endif
