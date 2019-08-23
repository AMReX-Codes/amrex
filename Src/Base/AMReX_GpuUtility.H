#ifndef AMREX_GPU_UTILITY_H_
#define AMREX_GPU_UTILITY_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <iostream>
#include <cmath>
#ifdef AMREX_USE_CUDA
#include <cuda.h>
#include <curand_kernel.h>
#endif

namespace amrex {
namespace Gpu {

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T LDG (Array4<T> const& a, int i, int j, int k) noexcept {
#ifdef __CUDA_ARCH__
        return __ldg(a.ptr(i,j,k,n));
#else
        return a(i,j,k);
#endif
    }

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T LDG (Array4<T> const& a, int i, int j, int k, int n) noexcept {
#ifdef __CUDA_ARCH__
        return __ldg(a.ptr(i,j,k,n));
#else
        return a(i,j,k,n);
#endif
    }

    inline bool isManaged (void const* p) noexcept {
#ifdef AMREX_USE_CUDA
        CUpointer_attribute attrib = CU_POINTER_ATTRIBUTE_IS_MANAGED;
        unsigned int is_managed = 0;
        void* data[] = { (void*)(&is_managed) };
        CUresult r = cuPointerGetAttributes(1, &attrib, data, (CUdeviceptr)p);
        return r == CUDA_SUCCESS && is_managed;
#else
        return false;
#endif
    }

    inline bool isDevicePtr (void const* p) noexcept {
#if defined(AMREX_USE_HIP)
        hipPointerAttribute_t attrib;
        // const_cast: earlier versions of hip have a bug
        hipError_t r = hipPointerGetAttributes(&attrib, const_cast<void*>(p));
        return r == hipSuccess && attrib.memoryType == hipMemoryTypeDevice;
#elif defined(AMREX_USE_CUDA)
        CUpointer_attribute attrib = CU_POINTER_ATTRIBUTE_MEMORY_TYPE;
        CUmemorytype mem_type = static_cast<CUmemorytype>(0);
        void* data[] = { (void*)(&mem_type) };
        CUresult r = cuPointerGetAttributes(1, &attrib, data, (CUdeviceptr)p);
        return r == CUDA_SUCCESS && mem_type == CU_MEMORYTYPE_DEVICE;
#else
        return false;
#endif
    }

    inline bool isPinnedPtr (void const* p) noexcept {
#if defined(AMREX_USE_HIP)
        hipPointerAttribute_t attrib;
        // const_cast: earlier versions of hip have a bug
        hipError_t r = hipPointerGetAttributes(&attrib, const_cast<void*>(p));
        return r == hipSuccess && attrib.memoryType == hipMemoryTypeHost;
#elif defined(AMREX_USE_CUDA)
        CUpointer_attribute attrib = CU_POINTER_ATTRIBUTE_MEMORY_TYPE;
        CUmemorytype mem_type = static_cast<CUmemorytype>(0);
        void* data[] = { (void*)(&mem_type) };
        CUresult r = cuPointerGetAttributes(1, &attrib, data, (CUdeviceptr)p);
        return r == CUDA_SUCCESS && mem_type == CU_MEMORYTYPE_HOST;
#else
        return false;
#endif
    }

    inline bool isGpuPtr (void const* p) noexcept {
#if defined(AMREX_USE_HIP)
        if (isManaged(p)) { // We might be using CUDA/NVCC
            return true;
        } else {
            hipPointerAttribute_t attrib;
            // const_cast: earlier versions of hip have a bug
            hipError_t r = hipPointerGetAttributes(&attrib, const_cast<void*>(p));
            return r == hipSuccess &&
                (attrib.memoryType == hipMemoryTypeHost   ||
                 attrib.memoryType == hipMemoryTypeDevice);
        }
#elif defined(AMREX_USE_CUDA)
        CUpointer_attribute attrib = CU_POINTER_ATTRIBUTE_MEMORY_TYPE;
        CUmemorytype mem_type = static_cast<CUmemorytype>(0);
        void* data[] = { (void*)(&mem_type) };
        CUresult r = cuPointerGetAttributes(1, &attrib, data, (CUdeviceptr)p);
        return r == CUDA_SUCCESS &&
            (mem_type == CU_MEMORYTYPE_HOST   ||
             mem_type == CU_MEMORYTYPE_DEVICE ||
             mem_type == CU_MEMORYTYPE_ARRAY  ||
             mem_type == CU_MEMORYTYPE_UNIFIED);
#else
        return false;
#endif
    }

    namespace Atomic {

#ifdef AMREX_USE_GPU
        namespace detail {

            AMREX_GPU_DEVICE AMREX_FORCE_INLINE
            float atomicMax(float* address, float val) noexcept
            {
                int* address_as_i = (int*) address;
                int old = *address_as_i, assumed;
                do {
                    assumed = old;
                    old = atomicCAS(address_as_i, assumed,
                                    __float_as_int(fmaxf(val, __int_as_float(assumed))));
                } while (assumed != old);
                return __int_as_float(old);
            }
            
            AMREX_GPU_DEVICE AMREX_FORCE_INLINE
            double atomicMax(double* address, double val) noexcept
            {
                unsigned long long int* address_as_ull = 
                    (unsigned long long int*) address;
                unsigned long long int old = *address_as_ull, assumed;
                do {
                    assumed = old;
                    old = atomicCAS(address_as_ull, assumed,
                                    __double_as_longlong(fmax(val, __longlong_as_double(assumed))));
                } while (assumed != old);
                return __longlong_as_double(old);
            }

            AMREX_GPU_DEVICE AMREX_FORCE_INLINE
            float atomicMin(float* address, float val) noexcept
            {
                int* address_as_i = (int*) address;
                int old = *address_as_i, assumed;
                do {
                    assumed = old;
                    old = atomicCAS(address_as_i, assumed,
                                    __float_as_int(fminf(val, __int_as_float(assumed))));
                } while (assumed != old);
                return __int_as_float(old);
            }

            AMREX_GPU_DEVICE AMREX_FORCE_INLINE
            double atomicMin(double* address, double val) noexcept
            {
                unsigned long long int* address_as_ull = 
                    (unsigned long long int*) address;
                unsigned long long int old = *address_as_ull, assumed;
                do {
                    assumed = old;
                    old = atomicCAS(address_as_ull, assumed,
                                    __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
                } while (assumed != old);
                return __longlong_as_double(old);
            }            
        } // namespace detail
#endif  // AMREX_USE_GPU

        template<class T>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Add (T* sum, T value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            atomicAdd(sum, value);
#else
            *sum += value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Add (long* sum, long value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            atomicAdd((unsigned long long*)sum, static_cast<unsigned long long>(value));
#else
            *sum += value;
#endif
        }
        
        template<class T>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Min (T* m, T value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            detail::atomicMin(m, value);
#else
            *m = (*m) < value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Min (int* m, int value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            atomicMin(m, value);
#else
            *m = (*m) < value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Min (unsigned int* m, unsigned int value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            atomicMin(m, value);
#else
            *m = (*m) < value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Min (unsigned long long int* m, unsigned long long int value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            atomicMin(m, value);
#else
            *m = (*m) < value ? (*m) : value;
#endif
        }
        
        template<class T>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Max (T* m, T value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            detail::atomicMax(m, value);
#else
            *m = (*m) > value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Max (int* m, int value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            atomicMax(m, value);
#else
            *m = (*m) > value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Max (unsigned int* m, unsigned int value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            atomicMax(m, value);
#else
            *m = (*m) > value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Max (unsigned long long int* m, unsigned long long int value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            atomicMax(m, value);
#else
            *m = (*m) > value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void LogicalOr (int* m, int value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            atomicOr(m, value);
#else
            *m = (*m) || value; 
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void LogicalAnd (int* m, int value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            atomicAnd(m, value ? ~0x0 : 0);
#else
            *m = (*m) && value; 
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        unsigned int Inc (unsigned int* m, unsigned int value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            return atomicInc(m, value);
#else
            return (*m)++;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        unsigned int Dec (unsigned int* m, unsigned int value) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            return atomicDec(m, value);
#else
            return (*m)--;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int Exch (int* address, int val) noexcept {
#ifdef AMREX_DEVICE_COMPILE
            return atomicExch(address, val);
#else
            int old = *address;
            *address = val;
            return old;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int CAS (int* address, int compare, int val) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            return atomicCAS(address, compare, val);
#else
            int old = *address;
            *address = (old == compare ? val : old);
            return old;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        unsigned int CAS (unsigned int* address,
                          unsigned int compare,
                          unsigned int val) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            return atomicCAS(address, compare, val);
#else
            unsigned int old = *address;
            *address = (old == compare ? val : old);
            return old;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        unsigned long long int CAS (unsigned long long int* address,
                                    unsigned long long int compare,
                                    unsigned long long int val) noexcept
        {
#ifdef AMREX_DEVICE_COMPILE
            return atomicCAS(address, compare, val);
#else
            unsigned long long int old = *address;
            *address = (old == compare ? val : old);
            return old;
#endif
        }
    } // namespace Atomic

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isnan (T m) noexcept
    {
#ifdef AMREX_DEVICE_COMPILE
        return ::isnan(m);
#else
        return std::isnan(m);
#endif
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isinf (T m) noexcept
    {
#ifdef AMREX_DEVICE_COMPILE
        return ::isinf(m);
#else
        return std::isinf(m);
#endif
    }

    struct StreamItInfo
    {
        bool device_sync;
        StreamItInfo () noexcept 
            : device_sync(true) {}
        StreamItInfo& DisableDeviceSync () noexcept {
            device_sync = false;
            return *this;
        }
    };

    class StreamIter
    {
    public:
        StreamIter (const int n, bool is_thread_safe=true) noexcept;
        StreamIter (const int n, const StreamItInfo& info, bool is_thread_safe=true) noexcept;

        ~StreamIter ();

        StreamIter (StreamIter const&) = delete;
        StreamIter (StreamIter &&) = delete;
        void operator= (StreamIter const&) = delete;
        void operator= (StreamIter &&) = delete;

        int operator() () const noexcept { return m_i; }

        bool isValid () const noexcept { return m_i < m_n; }

#if !defined(AMREX_USE_GPU)
        void operator++ () noexcept { ++m_i; }
#else
        void operator++ () noexcept;
#endif

    private:
        void init() noexcept;

        int m_n;
        int m_i;
        bool m_threadsafe;
        bool m_sync;
    };

} // namespace Gpu

#ifdef AMREX_USE_GPU
std::ostream& operator<< (std::ostream& os, const dim3& d);
#endif

using Gpu::isnan;
using Gpu::isinf;

namespace HostDevice { namespace Atomic {

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Add (T* sum, T value) noexcept
    {
#ifdef AMREX_DEVICE_COMPILE
        Gpu::Atomic::Add(sum,value);
#else
#ifdef _OPENMP
#pragma omp atomic update
#endif
        *sum += value;
#endif
    }

}}

} // namespace amrex

#endif
