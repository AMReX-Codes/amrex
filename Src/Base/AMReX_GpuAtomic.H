#ifndef AMREX_GPU_ATOMIC_H_
#define AMREX_GPU_ATOMIC_H_

#include <AMReX_GpuQualifiers.H>

namespace amrex {

namespace Gpu { namespace Atomic {

////////////////////////////////////////////////////////////////////////
//  Add
////////////////////////////////////////////////////////////////////////

    template<class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Add (T* sum, T value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        atomicAdd(sum, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic<T,as> a{sycl::multi_ptr<T,as>(sum)};
        sycl::atomic_fetch_add(a, value, mo);
#else
        *sum += value;
#endif
    }

#ifdef AMREX_USE_DPCPP

    namespace detail {
        template <typename R, typename I>
        AMREX_FORCE_INLINE
        void atomicAdd (R* address, R val) noexcept
        {
#if defined(__SYCL_DEVICE_ONLY__)
            constexpr auto mo = sycl::memory_order::relaxed;
            constexpr auto as = sycl::access::address_space::global_space;
            static_assert(sizeof(R) == sizeof(I), "sizeof R != sizeof I");
            I* add_as_I = reinterpret_cast<I*>(address);
            sycl::atomic<I,as> a{sycl::multi_ptr<I,as>(add_as_I)};
            I old_I = a.load(mo), new_I;
            do {
                R const old_R = *(reinterpret_cast<R const*>(&old_I));
                R const new_R = old_R + val;
                new_I = *(reinterpret_cast<I const*>(&new_R));
            } while (not a.compare_exchange_strong(old_I, new_I, mo));
#else
            *address += val;
#endif
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Add (double* sum, double value) noexcept
    {
        detail::atomicAdd<double, unsigned long long>(sum, value);
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Add (float* sum, float value) noexcept
    {
        detail::atomicAdd<float, int>(sum, value);
    }

#elif defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Add (long* sum, long value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        atomicAdd((unsigned long long*)sum, static_cast<unsigned long long>(value));
#else
        *sum += value;
#endif
    }

#endif

////////////////////////////////////////////////////////////////////////
//  Min
////////////////////////////////////////////////////////////////////////

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    namespace detail {
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        float atomicMin (float* address, float val) noexcept
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
        double atomicMin (double* address, double val) noexcept
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
    }
#endif

    template<class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Min (T* m, T value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        detail::atomicMin(m, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic<T,as> a{sycl::multi_ptr<T,as>(m)};
        sycl::atomic_fetch_min(a, value, mo);
#else
        *m = (*m) < value ? (*m) : value;
#endif
    }

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Min (int* m, int value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        atomicMin(m, value);
#else
        *m = (*m) < value ? (*m) : value;
#endif
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Min (unsigned int* m, unsigned int value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        atomicMin(m, value);
#else
        *m = (*m) < value ? (*m) : value;
#endif
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Min (unsigned long long int* m, unsigned long long int value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        atomicMin(m, value);
#else
        *m = (*m) < value ? (*m) : value;
#endif
    }

#endif

////////////////////////////////////////////////////////////////////////
//  Max
////////////////////////////////////////////////////////////////////////

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    namespace detail {
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        float atomicMax (float* address, float val) noexcept
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
        double atomicMax (double* address, double val) noexcept
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
    }
#endif

    template<class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Max (T* m, T value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        detail::atomicMax(m, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic<T,as> a{sycl::multi_ptr<T,as>(m)};
        sycl::atomic_fetch_max(a, value, mo);
#else
        *m = (*m) > value ? (*m) : value;
#endif
    }

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Max (int* m, int value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        atomicMax(m, value);
#else
        *m = (*m) > value ? (*m) : value;
#endif
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Max (unsigned int* m, unsigned int value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        atomicMax(m, value);
#else
        *m = (*m) > value ? (*m) : value;
#endif
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Max (unsigned long long int* m, unsigned long long int value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        atomicMax(m, value);
#else
        *m = (*m) > value ? (*m) : value;
#endif
    }

#endif

////////////////////////////////////////////////////////////////////////
//  LogicalOr
////////////////////////////////////////////////////////////////////////

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void LogicalOr (int* m, int value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        atomicOr(m, value);
#else
        *m = (*m) || value;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  LogicalAnd
////////////////////////////////////////////////////////////////////////

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void LogicalAnd (int* m, int value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        atomicAnd(m, value ? ~0x0 : 0);
#else
        *m = (*m) && value;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  Inc
////////////////////////////////////////////////////////////////////////

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    unsigned int Inc (unsigned int* m, unsigned int value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicInc(m, value);
#else
        return (*m)++;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  Dec
////////////////////////////////////////////////////////////////////////

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    unsigned int Dec (unsigned int* m, unsigned int value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicDec(m, value);
#else
        return (*m)--;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  Exch
////////////////////////////////////////////////////////////////////////

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T Exch (T* address, T val) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicExch(address, val);
#else
        T old = *address;
        *address = val;
        return old;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  CAS
////////////////////////////////////////////////////////////////////////

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T CAS (T* address, T compare, T val) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicCAS(address, compare, val);
#else
        T old = *address;
        *address = (old == compare ? val : old);
        return old;
#endif
    }
}}

namespace HostDevice { namespace Atomic {

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Add (T* sum, T value) noexcept
    {
#if AMREX_DEVICE_COMPILE()
        Gpu::Atomic::Add(sum,value);
#else
#ifdef _OPENMP
#pragma omp atomic update
#endif
        *sum += value;
#endif
    }

}}

}
#endif
