#ifndef AMREX_CUDAALLOCATORS_H_
#define AMREX_CUDAALLOCATORS_H_

#include <map>
#include <memory>
#include <limits>

#include <AMReX_Print.H>
#include <AMReX_Arena.H>
#include <AMReX_CudaDevice.H>

#ifdef AMREX_USE_CUDA
#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime.h>
#include <thrust/system/cuda/vector.h>
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
            cudaError_t error = cudaMallocManaged(&result, n*sizeof(T));
            if(error != cudaSuccess) {
                amrex::Print() << "Allocation failed in cudaMallocManaged \n";
                amrex::Abort();
            }

            const int device = Cuda::Device::deviceId();
            AMREX_GPU_SAFE_CALL(cudaMemAdvise(result, n*sizeof(T), 
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

    /**

       This is another managed memory allocator that inherits from the device allocator
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
	    cudaError_t error = cudaMallocManaged(&result, n*sizeof(T), cudaMemAttachGlobal);
  
	    if(error != cudaSuccess)
	    {
	        throw thrust::system_error(error, thrust::cuda_category(), 
					   "ThrustManagedAllocator::allocate(): cudaMallocManaged");
	    }
  
	    return thrust::device_pointer_cast(result);
	}
  
        inline void deallocate(pointer ptr, size_t)
        {
	    cudaError_t error = cudaFree(thrust::raw_pointer_cast(ptr));
  
	    if(error != cudaSuccess)
	    {
	        throw thrust::system_error(error, thrust::cuda_category(), 
					   "ThrustManagedAllocator::deallocate(): cudaFree");
	    }
	}
    };

    /**

       This allocator implements caching and is intended to be used internally by thrust when it allocates
       temporary storage.

       See:
           https://github.com/thrust/thrust/blob/master/examples/cuda/custom_temporary_allocation.cu

    **/
    class ThrustCachedAllocator
    {
    public:
        typedef char value_type;

        ThrustCachedAllocator() {}

        char *allocate(std::ptrdiff_t num_bytes)
        {
	  return (char*) The_Device_Arena()->alloc(num_bytes);
	}

        void deallocate(char *ptr, size_t n)
        {
	    The_Device_Arena()->free(ptr);
	}
      
    private:
        typedef std::multimap<std::ptrdiff_t, char*> free_blocks_type;
        typedef std::map<char *, std::ptrdiff_t>     allocated_blocks_type;
  
        free_blocks_type      free_blocks;
        allocated_blocks_type allocated_blocks; 
    };

    namespace Cuda
    {
        ThrustCachedAllocator& The_ThrustCachedAllocator ();
    }

// When not using CUDA, replace with standard allocator.
// Prevents need to wrap applications in lots of ifdefs.
#else

    template <typename T>
    using CudaManagedAllocator = std::allocator<T>;

#endif // AMREX_USE_CUDA

} // namespace amrex

#endif // AMREX_CUDAALLOCATORS_H_
