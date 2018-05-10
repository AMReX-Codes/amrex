#ifndef CUDAMANAGEDALLOCATOR_H_
#define CUDAMANAGEDALLOCATOR_H_

#include <limits>

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
        cudaError_t error = cudaMallocManaged(&result, cnt*sizeof(T), cudaMemAttachGlobal);
        if(error != cudaSuccess) {
            std::cout << "allocate failed in cudaMallocManaged" << std::endl;
        }        
        return result;
    }
    
    inline void deallocate(pointer p, size_type) {
        cudaError_t error = cudaFree(p);
        if(error != cudaSuccess) {
            std::cout << "deallocate failed in cudaFree" << std::endl;
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

#endif
