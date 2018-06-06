#ifndef CUDAALLOCATORS_H_
#define CUDAALLOCATORS_H_

#ifdef AMREX_USE_CUDA

#include <memory>
#include <iostream>
#include <limits>
#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime.h>
#include <AMReX_Print.H>
#include <AMReX.H>


namespace amrex {

template<typename T>
class CudaManagedAllocator {
public :

    typedef T value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    
    template<typename U>
    struct rebind {
        typedef CudaManagedAllocator<U> other;
    };
    
    inline explicit CudaManagedAllocator() {}
    inline         ~CudaManagedAllocator() {}
    inline explicit CudaManagedAllocator(CudaManagedAllocator const&) {}
    template<typename U>
    inline explicit CudaManagedAllocator(CudaManagedAllocator<U> const&) {}
    
    inline pointer       address(reference r)       { return &r; }
    inline const_pointer address(const_reference r) { return &r; }
    
    inline pointer allocate(size_type cnt,
                            typename std::allocator<void>::const_pointer = 0) {
        pointer result = nullptr;
        cudaError_t error = cudaMallocManaged(&result, cnt*sizeof(T));
        if(error != cudaSuccess) {
            amrex::Print() << "Allocation failed in cudaMallocManaged \n";
            amrex::Abort();
        }
        return result;
    }
    
    inline void deallocate(pointer p, size_type) {
        cudaError_t error = cudaFree(p);
        if(error != cudaSuccess) {
            amrex::Print() << "Deallocation failed in cudaFree \n";
            amrex::Abort();
        }
    }
    
    inline size_type max_size() const {
        return std::numeric_limits<size_type>::max() / sizeof(T);
    }
    
    inline void construct(pointer p, const T& t) { new(p) T(t); }
    inline void destroy(pointer p) { p->~T(); }
    
    inline bool operator==(CudaManagedAllocator const&)   { return true; }
    inline bool operator!=(CudaManagedAllocator const& a) { return !operator==(a); }
};

template<typename T>
class CudaHostAllocator {
public :

    typedef T value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    
    template<typename U>
    struct rebind {
        typedef CudaHostAllocator<U> other;
    };
    
    inline explicit CudaHostAllocator() {}
    inline         ~CudaHostAllocator() {}
    inline explicit CudaHostAllocator(CudaHostAllocator const&) {}
    template<typename U>
    inline explicit CudaHostAllocator(CudaHostAllocator<U> const&) {}
    
    inline pointer       address(reference r)       { return &r; }
    inline const_pointer address(const_reference r) { return &r; }
    
    inline pointer allocate(size_type cnt,
                            typename std::allocator<void>::const_pointer = 0) {
        pointer result = nullptr;
        cudaError_t error = cudaMallocHost(&result, cnt*sizeof(T));
        if(error != cudaSuccess) {
            amrex::Print() << "Allocation failed in cudaMallocHost \n";
            amrex::Abort();
        }
        return result;
    }
    
    inline void deallocate(pointer p, size_type) {
        cudaError_t error = cudaFreeHost(p);
        if(error != cudaSuccess) {
            amrex::Print() << "Deallocation failed in cudaFreeHost \n";
            amrex::Abort();
        }
    }
    
    inline size_type max_size() const {
        return std::numeric_limits<size_type>::max() / sizeof(T);
    }
    
    inline void construct(pointer p, const T& t) { new(p) T(t); }
    inline void destroy(pointer p) { p->~T(); }
    
    inline bool operator==(CudaHostAllocator const&)   { return true; }
    inline bool operator!=(CudaHostAllocator const& a) { return !operator==(a); }
};

}

#endif // CUDA
#endif // CUDAMANAGEDALLOCATOR_H_
