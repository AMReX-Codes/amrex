#ifndef AMREX_CUDAALLOCATORS_H_
#define AMREX_CUDAALLOCATORS_H_

#include <memory>
#include <limits>

#include <AMReX_Print.H>
#include <AMReX.H>

#ifdef AMREX_USE_CUDA
#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime.h>
#endif // AMREX_USE_CUDA

namespace amrex {

#ifdef AMREX_USE_CUDA

    template<typename T>
    class CudaManagedAllocator
    {
    public :

        using value_type = T;

        inline value_type* allocate(std::size_t n)
        {
            value_type* result = nullptr;
            cudaError_t error = cudaMallocManaged(&result, n*sizeof(T));
            if(error != cudaSuccess) {
                amrex::Print() << "Allocation failed in cudaMallocManaged \n";
                amrex::Abort();
            }

            const int device = Device::deviceId();
            CudaAPICheck(cudaMemAdvise(result, n*sizeof(T), 
                                       cudaMemAdviseSetPreferredLocation, device));

            return result;
        }
    
        inline void deallocate(value_type* p, std::size_t)
        {
            cudaError_t error = cudaFree(p);
            if(error != cudaSuccess) {
                amrex::Print() << "Deallocation failed in cudaFree \n";
                amrex::Abort();
            }
        }    
    };

    template <class T, class U>
    bool
    operator==(CudaManagedAllocator<T> const&, CudaManagedAllocator<U> const&) noexcept
    {
        return true;
    }
    
    template <class T, class U>
    bool
    operator!=(CudaManagedAllocator<T> const& x, CudaManagedAllocator<U> const& y) noexcept
    {
        return !(x == y);
    }
    
#endif // AMREX_USE_CUDA

} // namespace amrex

#endif // AMREX_CUDAALLOCATORS_H_
