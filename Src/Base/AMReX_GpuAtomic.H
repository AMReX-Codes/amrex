#ifndef AMREX_GPU_ATOMIC_H_
#define AMREX_GPU_ATOMIC_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Functional.H>

namespace amrex {

namespace Gpu { namespace Atomic {

// For Add, Min and Max, we support int, unsigned int, long, unsigned long long, float and double.
// For LogicalOr and LogicalAnd, the data type is int.
// For Inc and Dec, the data type is unsigned int.
// For Exch and CAS, the data type is generic.
// All these functions are non-atomic in host code!!!
// If one needs them to be atomic in host code, use HostDevice::Atomic::*.  Currently only
// HostDevice::Atomic is supported.  We could certainly add more.

namespace detail {

#ifdef AMREX_USE_DPCPP

    template <typename R, typename I, typename F>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    R atomic_op (R* const address, R const val, F const f) noexcept
    {
#if defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        static_assert(sizeof(R) == sizeof(I), "sizeof R != sizeof I");
        I* const add_as_I = reinterpret_cast<I*>(address);
        sycl::atomic<I,as> a{sycl::multi_ptr<I,as>(add_as_I)};
        I old_I = a.load(mo), new_I;
        do {
            R const new_R = f(*(reinterpret_cast<R const*>(&old_I)), val);
            new_I = *(reinterpret_cast<I const*>(&new_R));
        } while (not a.compare_exchange_strong(old_I, new_I, mo));
        return *(reinterpret_cast<R const*>(&old_I));
#else
        R old = *address;
        *address = f(*address, val);
        return old;
#endif
    }

#elif defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    template <typename R, typename I, typename F>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    R atomic_op (R* const address, R const val, F const f) noexcept
    {
        static_assert(sizeof(R) == sizeof(I), "sizeof R != sizeof I");
        I* const add_as_I = reinterpret_cast<I*>(address);
        I old_I = *add_as_I, assumed_I;
        do {
            assumed_I = old_I;
            R const new_R = f(*(reinterpret_cast<R const*>(&assumed_I)), val);
            old_I = atomicCAS(add_as_I, assumed_I, *(reinterpret_cast<I const*>(&new_R)));
        } while (assumed_I != old_I);
        return *(reinterpret_cast<R const*>(&old_I));
    }

#endif

}

////////////////////////////////////////////////////////////////////////
//  Add
////////////////////////////////////////////////////////////////////////

#ifdef AMREX_USE_GPU

    template<class T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T Add_device (T* const sum, T const value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicAdd(sum, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic<T,as> a{sycl::multi_ptr<T,as>(sum)};
        return sycl::atomic_fetch_add(a, value, mo);
#else
        amrex::ignore_unused(sum, value);
        return T(); // should never get here, but have to return something
#endif
    }

#ifdef AMREX_USE_DPCPP

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    float Add_device (float* const sum, float const value) noexcept
    {
        return detail::atomic_op<float, int>(sum, value, amrex::Plus<float>());
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    double Add_device (double* const sum, double const value) noexcept
    {
        return detail::atomic_op<double, unsigned long long>(sum, value, amrex::Plus<double>());
    }

#elif defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Long Add_device (Long* const sum, Long const value) noexcept
    {
        return detail::atomic_op<Long, unsigned long long>(sum, value, amrex::Plus<Long>());
    }

#endif

#endif

    template<class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T Add (T* const sum, T const value) noexcept
    {
#if AMREX_DEVICE_COMPILE
        return Add_device(sum, value);
#else
        auto old = *sum;
        *sum += value;
        return old;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  Min
////////////////////////////////////////////////////////////////////////

#ifdef AMREX_USE_GPU

    template<class T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T Min_device (T* const m, T const value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicMin(m, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic<T,as> a{sycl::multi_ptr<T,as>(m)};
        return sycl::atomic_fetch_min(a, value, mo);
#else
        amrex::ignore_unused(m,value);
        return T(); // should never get here, but have to return something
#endif
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    float Min_device (float* const m, float const value) noexcept
    {
        return detail::atomic_op<float,int>(m,value,amrex::Less<float>());
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    double Min_device (double* const m, double const value) noexcept
    {
        return detail::atomic_op<double, unsigned long long int>(m,value,amrex::Less<double>());
    }

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Long Min_device (Long* const m, Long const value) noexcept
    {
        return detail::atomic_op<Long, unsigned long long int>(m,value,amrex::Less<Long>());
    }

#endif

#endif

    template<class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T Min (T* const m, T const value) noexcept
    {
#if AMREX_DEVICE_COMPILE
        return Min_device(m, value);
#else
        auto const old = *m;
        *m = (*m) < value ? (*m) : value;
        return old;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  Max
////////////////////////////////////////////////////////////////////////

#ifdef AMREX_USE_GPU

    template<class T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T Max_device (T* const m, T const value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicMax(m, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic<T,as> a{sycl::multi_ptr<T,as>(m)};
        return sycl::atomic_fetch_max(a, value, mo);
#else
        amrex::ignore_unused(m,value);
        return T(); // should never get here, but have to return something
#endif
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    float Max_device (float* const m, float const value) noexcept
    {
        return detail::atomic_op<float,int>(m,value,amrex::Greater<float>());
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    double Max_device (double* const m, double const value) noexcept
    {
        return detail::atomic_op<double, unsigned long long int>(m,value,amrex::Greater<double>());
    }

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Long Max_device (Long* const m, Long const value) noexcept
    {
        return detail::atomic_op<Long, unsigned long long int>(m,value,amrex::Greater<Long>());
    }

#endif

#endif

    template<class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T Max (T* const m, T const value) noexcept
    {
#if AMREX_DEVICE_COMPILE
        return Max_device(m, value);
#else
        auto const old = *m;
        *m = (*m) > value ? (*m) : value;
        return old;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  LogicalOr
////////////////////////////////////////////////////////////////////////

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int LogicalOr (int* const m, int const value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicOr(m, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic<int,as> a{sycl::multi_ptr<int,as>(m)};
        return sycl::atomic_fetch_or(a, value, mo);
#else
        int const old = *m;
        *m = (*m) || value;
        return old;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  LogicalAnd
////////////////////////////////////////////////////////////////////////

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int LogicalAnd (int* const m, int const value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicAnd(m, value ? ~0x0 : 0);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic<int,as> a{sycl::multi_ptr<int,as>(m)};
        return sycl::atomic_fetch_and(a, value ? ~0x0 : 0, mo);
#else
        int const old = *m;
        *m = (*m) && value;
        return old;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  Inc
////////////////////////////////////////////////////////////////////////

#ifdef AMREX_USE_DPCPP
    template <sycl::access::address_space AS = sycl::access::address_space::global_space>
    AMREX_FORCE_INLINE
    unsigned int Inc (unsigned int* const m, unsigned int const value) noexcept
    {
#if defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        sycl::atomic<unsigned int,AS> a{sycl::multi_ptr<unsigned int,AS>(m)};
        unsigned int oldi = a.load(mo), newi;
        do {
            newi = (oldi >= value) ? 0u : (oldi+1u);
        } while (not a.compare_exchange_strong(oldi, newi, mo));
        return oldi;
#else
        auto const old = *m;
        *m = (old >= value) ? 0u : (old+1u);
        return old;
#endif
    }
#else
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    unsigned int Inc (unsigned int* const m, unsigned int const value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicInc(m, value);
#else
        auto const old = *m;
        *m = (old >= value) ? 0u : (old+1u);
        return old;
#endif
    }
#endif

////////////////////////////////////////////////////////////////////////
//  Dec
////////////////////////////////////////////////////////////////////////

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    unsigned int Dec (unsigned int* const m, unsigned int const value) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicDec(m, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic<unsigned int,as> a{sycl::multi_ptr<unsigned int,as>(m)};
        unsigned int oldi = a.load(mo), newi;
        do {
            newi = ((oldi == 0u) || (oldi > value)) ? value : (oldi-1u);
        } while (not a.compare_exchange_strong(oldi, newi, mo));
        return oldi;
#else
        auto const old = *m;
        *m = ((old == 0u) || (old > value)) ? value : (old-1u);
        return old;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  Exch
////////////////////////////////////////////////////////////////////////

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T Exch (T* const address, T const val) noexcept
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicExch(address, val);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic<T,as> a{sycl::multi_ptr<T,as>(address)};
        return sycl::atomic_exchange(a, val, mo);
#else
        auto const old = *address;
        *address = val;
        return old;
#endif
    }

////////////////////////////////////////////////////////////////////////
//  CAS
////////////////////////////////////////////////////////////////////////

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T CAS (T* const address, T compare, T const val) noexcept
    {           // cannot be T const compare because of compare_exchange_strong
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return atomicCAS(address, compare, val);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic<T,as> a{sycl::multi_ptr<T,as>(address)};
        a.compare_exchange_strong(compare, val, mo);
        return compare;
#else
        auto const old = *address;
        *address = (old == compare ? val : old);
        return old;
#endif
    }
}}

namespace HostDevice { namespace Atomic {

    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Add (T* const sum, T const value) noexcept
    {
#if AMREX_DEVICE_COMPILE
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
