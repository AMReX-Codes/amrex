#ifndef AMREX_GPU_ATOMIC_H_
#define AMREX_GPU_ATOMIC_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Functional.H>
#include <AMReX_INT.H>

#include <utility>

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
        constexpr auto ms = sycl::memory_scope::device;
        constexpr auto as = sycl::access::address_space::global_space;
        static_assert(sizeof(R) == sizeof(I), "sizeof R != sizeof I");
        I* const add_as_I = reinterpret_cast<I*>(address);
        sycl::atomic_ref<I,mo,ms,as> a{*add_as_I};
        I old_I = a.load(), new_I;
        do {
            R const new_R = f(*(reinterpret_cast<R const*>(&old_I)), val);
            new_I = *(reinterpret_cast<I const*>(&new_R));
        } while (! a.compare_exchange_strong(old_I, new_I));
        return *(reinterpret_cast<R const*>(&old_I));
#else
        R old = *address;
        *address = f(*address, val);
        return old;
#endif
    }

    template <typename R, typename I, typename Op, typename Cond>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool atomic_op_if (R* const address, R const val, Op&& op, Cond&& cond) noexcept
    {
#if defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto ms = sycl::memory_scope::device;
        constexpr auto as = sycl::access::address_space::global_space;
        static_assert(sizeof(R) == sizeof(I), "sizeof R != sizeof I");
        I* const add_as_I = reinterpret_cast<I*>(address);
        sycl::atomic_ref<I,mo,ms,as> a{*add_as_I};
        I old_I = a.load(), new_I;
        bool test_success;
        do {
            R const tmp = op(*(reinterpret_cast<R const*>(&old_I)), val);
            new_I = *(reinterpret_cast<I const*>(&tmp));
            test_success = cond(tmp);
        } while (test_success && ! a.compare_exchange_strong(old_I, new_I));
        return test_success;
#else
        R old = *address;
        R const tmp = op(*(reinterpret_cast<R const*>(&old)), val);
        if (cond(tmp)) {
            *address = tmp;
            return true;
        } else {
            return false;
        }
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

    template <typename R, typename I, typename Op, typename Cond>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool atomic_op_if (R* const address, R const val, Op&& op, Cond&& cond) noexcept
    {
        static_assert(sizeof(R) == sizeof(I), "sizeof R != sizeof I");
        I* const add_as_I = reinterpret_cast<I*>(address);
        I old_I = *add_as_I, assumed_I;
        bool test_success;
        do {
            assumed_I = old_I;
            R const new_R = op(*(reinterpret_cast<R const*>(&assumed_I)), val);
            test_success = cond(new_R);
            if (test_success) {
                old_I = atomicCAS(add_as_I, assumed_I, *(reinterpret_cast<I const*>(&new_R)));
            }
        } while (test_success && assumed_I != old_I);
        return test_success;
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
#if defined(__CUDA_ARCH__) && defined(AMREX_USE_CUDA) || defined(__HIP_DEVICE_COMPILE__)
        return atomicAdd(sum, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto ms = sycl::memory_scope::device;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic_ref<T,mo,ms,as> a{*sum};
        return a.fetch_add(value);
#else
        amrex::ignore_unused(sum, value);
        return T(); // should never get here, but have to return something
#endif
    }

#if defined(AMREX_USE_HIP) && defined(__gfx90a__)
    // https://github.com/ROCm-Developer-Tools/hipamd/blob/rocm-4.5.x/include/hip/amd_detail/amd_hip_unsafe_atomics.h
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    float Add_device (float* const sum, float const value) noexcept
    {
        return unsafeAtomicAdd(sum, value);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    double Add_device (double* const sum, double const value) noexcept
    {
        return unsafeAtomicAdd(sum, value);
    }
#endif

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    // https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomicadd
    // https://rocmdocs.amd.com/en/latest/Programming_Guides/Kernel_language.html?#atomic-functions
    // CUDA and HIP support int, unsigned int, and unsigned long long.

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Long Add_device (Long* const sum, Long const value) noexcept
    {
        return detail::atomic_op<Long, unsigned long long>(sum, value, amrex::Plus<Long>());
    }

#endif

#if defined(AMREX_USE_CUDA) && (__CUDA_ARCH__ < 600)

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    double Add_device (double* const sum, double const value) noexcept
    {
        return detail::atomic_op<double, unsigned long long>(sum, value, amrex::Plus<double>());
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
//  If
////////////////////////////////////////////////////////////////////////

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP) || defined(AMREX_USE_DPCPP)
    template <typename T, typename Op, typename Cond,
              std::enable_if_t<sizeof(T) == sizeof(unsigned int), int> foo = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool If_device (T* const sum, T const value, Op&& op, Cond&& cond) noexcept
    {
        return detail::atomic_op_if<T, unsigned int>(sum, value,
                   std::forward<Op>(op), std::forward<Cond>(cond));
    }

    template <typename T, typename Op, typename Cond,
              std::enable_if_t<sizeof(T) == sizeof(unsigned long long), int> foo = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    bool If_device (T* const sum, T const value, Op&& op, Cond&& cond) noexcept
    {
        return detail::atomic_op_if<T, unsigned long long>(sum, value,
                   std::forward<Op>(op), std::forward<Cond>(cond));
    }
#endif

/**
*  \brief Conditionally perform an atomic operation.
*
*  Atomically updates the result at "add" with "value" using "op",
*  but only if "cond" is true.
*
* \tparam T the type pointed to by add
* \tparam Op callable that takes two T argument and combines them
* \tparam Cond callable that takes a "T" a returns whether to do the update
*
* \param add address to atomically update
* \param value value to combine
* \param op callable specifying the operation to use to combine *add and value
* \param cond callable specifying the condition to test on first.
*        The value passed in to the cond function is the would-be combined value
**/
    template<class T, class Op, class Cond>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool If (T* const add, T const value, Op&& op, Cond&& cond) noexcept
    {
#if AMREX_DEVICE_COMPILE
        return If_device(add, value, std::forward<Op>(op), std::forward<Cond>(cond));
#else
        T old = *add;
        T const tmp = op(old, value);
        if (cond(tmp)) {
            *add = tmp;
            return true;
        } else {
            return false;
        }
#endif
    }

////////////////////////////////////////////////////////////////////////
//  AddNoRet
////////////////////////////////////////////////////////////////////////

    template<class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void AddNoRet (T* const sum, T const value) noexcept
    {
#if AMREX_DEVICE_COMPILE
        Add_device(sum, value);
#else
        *sum += value;
#endif
    }

#if defined(AMREX_USE_HIP) && defined(HIP_VERSION_MAJOR) && (HIP_VERSION_MAJOR < 5)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void AddNoRet (float* const sum, float const value) noexcept
    {
#if AMREX_DEVICE_COMPILE
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
        atomicAddNoRet(sum, value);
#pragma clang diagnostic pop
#else
        *sum += value;
#endif
    }
#endif

////////////////////////////////////////////////////////////////////////
//  Min
////////////////////////////////////////////////////////////////////////

#ifdef AMREX_USE_GPU

    template<class T>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    T Min_device (T* const m, T const value) noexcept
    {
#if defined(__CUDA_ARCH__) && defined(AMREX_USE_CUDA) || defined(__HIP_DEVICE_COMPILE__)
        return atomicMin(m, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto ms = sycl::memory_scope::device;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic_ref<T,mo,ms,as> a{*m};
        return a.fetch_min(value);
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
#if defined(__CUDA_ARCH__) && defined(AMREX_USE_CUDA) || defined(__HIP_DEVICE_COMPILE__)
        return atomicMax(m, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto ms = sycl::memory_scope::device;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic_ref<T,mo,ms,as> a{*m};
        return a.fetch_max(value);
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
#if defined(__CUDA_ARCH__) && defined(AMREX_USE_CUDA) || defined(__HIP_DEVICE_COMPILE__)
        return atomicOr(m, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto ms = sycl::memory_scope::device;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic_ref<int,mo,ms,as> a{*m};
        return a.fetch_or(value);
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
#if defined(__CUDA_ARCH__) && defined(AMREX_USE_CUDA) || defined(__HIP_DEVICE_COMPILE__)
        return atomicAnd(m, value ? ~0x0 : 0);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto ms = sycl::memory_scope::device;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic_ref<int,mo,ms,as> a{*m};
        return a.fetch_and(value ? ~0x0 : 0);
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
        constexpr auto ms = sycl::memory_scope::device;
        sycl::atomic_ref<unsigned int,mo,ms,AS> a{*m};
        unsigned int oldi = a.load(), newi;
        do {
            newi = (oldi >= value) ? 0u : (oldi+1u);
        } while (! a.compare_exchange_strong(oldi, newi));
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
#if defined(__CUDA_ARCH__) && defined(AMREX_USE_CUDA) || defined(__HIP_DEVICE_COMPILE__)
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
#if defined(__CUDA_ARCH__) && defined(AMREX_USE_CUDA) || defined(__HIP_DEVICE_COMPILE__)
        return atomicDec(m, value);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto ms = sycl::memory_scope::device;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic_ref<unsigned int,mo,ms,as> a{*m};
        unsigned int oldi = a.load(), newi;
        do {
            newi = ((oldi == 0u) || (oldi > value)) ? value : (oldi-1u);
        } while (! a.compare_exchange_strong(oldi, newi));
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
#if defined(__CUDA_ARCH__) && defined(AMREX_USE_CUDA) || defined(__HIP_DEVICE_COMPILE__)
        return atomicExch(address, val);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto ms = sycl::memory_scope::device;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic_ref<T,mo,ms,as> a{*address};
        return a.exchange(val);
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
#if defined(__CUDA_ARCH__) && defined(AMREX_USE_CUDA) || defined(__HIP_DEVICE_COMPILE__)
        return atomicCAS(address, compare, val);
#elif defined(__SYCL_DEVICE_ONLY__)
        constexpr auto mo = sycl::memory_order::relaxed;
        constexpr auto ms = sycl::memory_scope::device;
        constexpr auto as = sycl::access::address_space::global_space;
        sycl::atomic_ref<T,mo,ms,as> a{*address};
        a.compare_exchange_strong(compare, val);
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
        Gpu::Atomic::AddNoRet(sum,value);
#else
#ifdef AMREX_USE_OMP
#pragma omp atomic update
#endif
        *sum += value;
#endif
    }

}}

#ifdef AMREX_USE_GPU
// functors
namespace Gpu {
    template <typename T>
    struct AtomicAdd
    {
        AMREX_GPU_DEVICE void operator() (T* const dest, T const source) noexcept {
            Gpu::Atomic::AddNoRet(dest, source);
        }
    };

    template <typename T>
    struct AtomicMin
    {
        AMREX_GPU_DEVICE void operator() (T* const dest, T const source) noexcept {
            Gpu::Atomic::Min(dest, source);
        }
    };

    template <typename T>
    struct AtomicMax
    {
        AMREX_GPU_DEVICE void operator() (T* const dest, T const source) noexcept {
            Gpu::Atomic::Max(dest, source);
        }
    };

    template <typename T>
    struct AtomicLogicalAnd
    {
        AMREX_GPU_DEVICE void operator() (T* const dest, T const source) noexcept {
            Gpu::Atomic::LogicalAnd(dest, source);
        }
    };

    template <typename T>
    struct AtomicLogicalOr
    {
        AMREX_GPU_DEVICE void operator() (T* const dest, T const source) noexcept {
            Gpu::Atomic::LogicalOr(dest, source);
        }
    };
}
#endif

}
#endif
