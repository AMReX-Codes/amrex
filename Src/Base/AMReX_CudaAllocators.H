#ifndef AMREX_CUDAALLOCATORS_H_
#define AMREX_CUDAALLOCATORS_H_

#include <map>
#include <memory>
#include <limits>

#include <AMReX_Print.H>
#include <AMReX_Arena.H>
#include <AMReX_GpuDevice.H>

#ifdef AMREX_USE_CUDA
#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime.h>
#include <thrust/system/cuda/vector.h>
#include <thrust/device_malloc_allocator.h>
#endif // AMREX_USE_CUDA

namespace amrex {

#ifdef AMREX_USE_CUDA
  
    /**

       This is a custom memory allocator that wraps around cudaMallocManaged.
       Intended for use with std::vector.

    **/
    template<typename T>
    class CudaManagedAllocator
    {
    public :

        using value_type = T;

        inline value_type* allocate(std::size_t n)
        {
            value_type* result = nullptr;
            AMREX_CUDA_SAFE_CALL(cudaMallocManaged(&result, n*sizeof(T)));

            const int device = Gpu::Device::deviceId();
            AMREX_CUDA_SAFE_CALL(cudaMemAdvise(result, n*sizeof(T), 
                                              cudaMemAdviseSetPreferredLocation, device));

            return result;
        }
    
        inline void deallocate(value_type* p, std::size_t)
        {
            AMREX_CUDA_SAFE_CALL(cudaFree(p));
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

    /**

       This is a custom memory allocator that wraps around cudaMallocManaged.
       Intended for use with std::vector.

    **/
    template<typename T>
    class CudaPinnedAllocator
    {
    public :

        using value_type = T;

        inline value_type* allocate(std::size_t n)
        {
            value_type* result = nullptr;
            AMREX_CUDA_SAFE_CALL(cudaMallocHost(&result, n*sizeof(T)));

            return result;
        }
    
        inline void deallocate(value_type* p, std::size_t)
        {
            AMREX_CUDA_SAFE_CALL(cudaFreeHost(p));
        }    
    };

    template <class T, class U>
    bool
    operator==(CudaPinnedAllocator<T> const&, CudaPinnedAllocator<U> const&) noexcept
    {
        return true;
    }
    
    template <class T, class U>
    bool
    operator!=(CudaPinnedAllocator<T> const& x, CudaPinnedAllocator<U> const& y) noexcept
    {
        return !(x == y);
    }

    /**

       This is a managed memory allocator that inherits from the device allocator
       in thrust. Intended for use with thrust::device_vector.

    **/
    template<class T>
    class ThrustManagedAllocator : public thrust::device_malloc_allocator<T>
    {
    public:
        using value_type = T;
    
        typedef thrust::device_ptr<T>  pointer;
        inline pointer allocate(size_t n)
        {
	    value_type* result = nullptr;
            result = (value_type*) The_Arena()->alloc(n * sizeof(T));
	    return thrust::device_pointer_cast(result);
	}
  
        inline void deallocate(pointer ptr, size_t)
        {
            The_Arena()->free(thrust::raw_pointer_cast(ptr));
	}
    };

    /**

       This is a memory allocator that uses the Pinned memory Pool from AMReX.
       It is intended for use with thrust::device_vector.

    **/
    template<class T>
    class ThrustPinnedAllocator : public thrust::device_malloc_allocator<T>
    {
    public:
        using value_type = T;
    
        typedef thrust::device_ptr<T>  pointer;
        inline pointer allocate(size_t n)
        {
	    value_type* result = nullptr;
            result = (value_type*) The_Pinned_Arena()->alloc(n * sizeof(T));  
	    return thrust::device_pointer_cast(result);
	}
  
        inline void deallocate(pointer ptr, size_t)
        {
            The_Pinned_Arena()->free(thrust::raw_pointer_cast(ptr));
	}
    };

    /**

       This is a memory allocator that uses the Device memory Pool from AMReX.
       It is intended for use with thrust::device_vector.

    **/
    template<class T>
    class ThrustDeviceAllocator : public thrust::device_malloc_allocator<T>
    {
    public:
        using value_type = T;
    
        typedef thrust::device_ptr<T>  pointer;
        inline pointer allocate(size_t n)
        {
	    value_type* result = nullptr;
            result = (value_type*) The_Device_Arena()->alloc(n * sizeof(T));  
	    return thrust::device_pointer_cast(result);
	}
  
        inline void deallocate(pointer ptr, size_t)
        {
            The_Device_Arena()->free(thrust::raw_pointer_cast(ptr));
	}
    };

    /**

       This is a polymorphic memory allocator that uses either the 
       Pinned or Device memory pool from AMReX, depending on the value of
       amrex.use_gpu_aware_mpi. It is intended for use with thrust::device_vector.

    **/
    template<class T>
    class ThrustPolymorphicAllocator : public thrust::device_malloc_allocator<T>
    {
    public:
        using value_type = T;
    
        typedef thrust::device_ptr<T>  pointer;
        inline pointer allocate(size_t n)
        {
	    value_type* result = nullptr;
            if (ParallelDescriptor::UseGpuAwareMpi()) 
            {
                result = (value_type*) The_Device_Arena()->alloc(n * sizeof(T));
            }
            else
            {
                result = (value_type*) The_Pinned_Arena()->alloc(n * sizeof(T));
            }
	    return thrust::device_pointer_cast(result);
	}
  
        inline void deallocate(pointer ptr, size_t)
        {
            if (ParallelDescriptor::UseGpuAwareMpi()) 
            {
                The_Device_Arena()->free(thrust::raw_pointer_cast(ptr));
            }
            else
            {
                The_Pinned_Arena()->free(thrust::raw_pointer_cast(ptr));
            }
	}
    };

    namespace Cuda
    {
        ThrustManagedAllocator<char>& The_ThrustCachedAllocator ();
    }

// When not using CUDA, replace with standard allocator.
// Prevents need to wrap applications in lots of ifdefs.
#else

    template <typename T>
    using CudaManagedAllocator = std::allocator<T>;

#endif // AMREX_USE_CUDA

} // namespace amrex

#endif // AMREX_CUDAALLOCATORS_H_
