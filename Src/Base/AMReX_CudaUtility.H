#ifndef AMREX_CUDA_UTILITY_H_
#define AMREX_CUDA_UTILITY_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_CudaDevice.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <iostream>
#include <cmath>
#ifdef AMREX_USE_CUDA
#include <cuda.h>
#include <curand_kernel.h>
#endif

namespace amrex {
namespace Cuda {

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
#ifdef AMREX_USE_CUDA
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
#ifdef AMREX_USE_CUDA
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
#ifdef AMREX_USE_CUDA
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

#ifdef AMREX_USE_CUDA
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
#endif  // AMREX_USE_CUDA

        template<class T>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Add (T* sum, T value) noexcept
        {
#if defined(__CUDA_ARCH__)
            atomicAdd(sum, value);
#else
            *sum += value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Add (long* sum, long value) noexcept
        {
#if defined(__CUDA_ARCH__)
            atomicAdd((unsigned long long*)sum, static_cast<unsigned long long>(value));
#else
            *sum += value;
#endif
        }
        
        template<class T>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Min (T* m, T value) noexcept
        {
#if defined(__CUDA_ARCH__)
            detail::atomicMin(m, value);
#else
            *m = (*m) < value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Min (int* m, int value) noexcept
        {
#if defined(__CUDA_ARCH__)
            atomicMin(m, value);
#else
            *m = (*m) < value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Min (unsigned int* m, unsigned int value) noexcept
        {
#if defined(__CUDA_ARCH__)
            atomicMin(m, value);
#else
            *m = (*m) < value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Min (unsigned long long int* m, unsigned long long int value) noexcept
        {
#if defined(__CUDA_ARCH__)
            atomicMin(m, value);
#else
            *m = (*m) < value ? (*m) : value;
#endif
        }
        
        template<class T>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Max (T* m, T value) noexcept
        {
#if defined(__CUDA_ARCH__)
            detail::atomicMax(m, value);
#else
            *m = (*m) > value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Max (int* m, int value) noexcept
        {
#if defined(__CUDA_ARCH__)
            atomicMax(m, value);
#else
            *m = (*m) > value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Max (unsigned int* m, unsigned int value) noexcept
        {
#if defined(__CUDA_ARCH__)
            atomicMax(m, value);
#else
            *m = (*m) > value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Max (unsigned long long int* m, unsigned long long int value) noexcept
        {
#if defined(__CUDA_ARCH__)
            atomicMax(m, value);
#else
            *m = (*m) > value ? (*m) : value;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void Or (int* m, int value) noexcept
        {
#if defined(__CUDA_ARCH__)
            atomicOr(m, value);
#else
            *m = (*m) || value; 
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void And (int* m, int value) noexcept
        {
#if defined(__CUDA_ARCH__)
            atomicAnd(m, value);
#else
            *m = (*m) && value; 
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        unsigned int Inc (unsigned int* m, unsigned int value) noexcept
        {
#if defined(__CUDA_ARCH__)
            return atomicInc(m, value);
#else
            return (*m)++;
#endif
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        unsigned int Dec (unsigned int* m, unsigned int value) noexcept
        {
#if defined(__CUDA_ARCH__)
            return atomicDec(m, value);
#else
            return (*m)--;
#endif
        }

    } // namespace Atomic

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isnan(T m) noexcept
    {
#if defined(__CUDA_ARCH__)
        return isnan(m);
#else
        return std::isnan(m);
#endif
    }

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isinf(T m) noexcept
    {
#ifdef __CUDA_ARCH__
        return isinf(m);
#else
        return std::isinf(m);
#endif
    }

    class StreamIter
    {
    public:
        StreamIter (const int n, bool is_thread_safe=true) noexcept;
        ~StreamIter ();

        StreamIter (StreamIter const&) = delete;
        StreamIter (StreamIter &&) = delete;
        void operator= (StreamIter const&) = delete;
        void operator= (StreamIter &&) = delete;

        int operator() () const noexcept { return m_i; }

        bool isValid () const noexcept { return m_i < m_n; }

#if !defined(AMREX_USE_CUDA)
        void operator++ () noexcept { ++m_i; }
#else
        void operator++ () noexcept;
#endif

    private:
        int m_n;
        int m_i;
        bool m_threadsafe;
    };

} // namespace Cuda

#ifdef AMREX_USE_CUDA
std::ostream& operator<< (std::ostream& os, const dim3& d);
#endif

using Cuda::isnan;
using Cuda::isinf;

} // namespace amrex

#endif
