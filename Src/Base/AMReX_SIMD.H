#ifndef AMREX_SIMD_H_
#define AMREX_SIMD_H_

#include <AMReX_Config.H>

#include <AMReX.H>  // amrex::ignore_unused
#include <AMReX_REAL.H>

#ifdef AMREX_USE_SIMD
// TODO make SIMD provider configurable: VIR (C++17 TS2) or C++26 (later)
#   include <vir/simd.h>  // includes SIMD TS2 header <experimental/simd>
#   if __cplusplus >= 202002L
#       include <vir/simd_cvt.h>
#   endif
#endif

#include <cstdint>
#include <type_traits>


namespace amrex::simd
{
    // TODO make SIMD provider configurable: VIR (C++17 TS2) or C++26 (later)
    // for https://en.cppreference.com/w/cpp/experimental/simd/simd_cast.html
    namespace stdx {
#ifdef AMREX_USE_SIMD
        using namespace vir::stdx;
#   if __cplusplus >= 202002L
        using vir::cvt;
#   endif
#else
        // fallback implementations for functions that are commonly used in portable code paths

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool any_of (bool const v) { return v; }
#endif
    }

    // TODO: move to AMReX_REAL.H?

#ifdef AMREX_USE_SIMD
    // TODO: not sure why std::experimental::simd_abi::native<T> does not work, so we use this long version
    constexpr auto native_simd_size_real = stdx::native_simd<amrex::Real>::size();
    constexpr auto native_simd_size_particlereal = stdx::native_simd<amrex::ParticleReal>::size();

    // Note: to make use of not only vector registers but also ILP, user might want to use * 2 or more of the native size
    //       for selected compute kernels.
    // TODO Check if a default with * 2 or similar is sensible.
    template<int SIMD_WIDTH = native_simd_size_real>
    using SIMDReal = stdx::fixed_size_simd<amrex::Real, SIMD_WIDTH>;

    template<int SIMD_WIDTH = native_simd_size_particlereal>
    using SIMDParticleReal = stdx::fixed_size_simd<amrex::ParticleReal, SIMD_WIDTH>;

    // Type that has the same number of int SIMD elements as the SIMDParticleReal type
    template<typename T_ParticleReal = SIMDParticleReal<>>
    using SIMDInt = stdx::rebind_simd_t<int, T_ParticleReal>;

    // Type that has the same number of IdCpu SIMD elements as the SIMDParticleReal type
    template<typename T_ParticleReal = SIMDParticleReal<>>
    using SIMDIdCpu = stdx::rebind_simd_t<std::uint64_t, T_ParticleReal>;
#else
    constexpr auto native_simd_size_real = 1;
    constexpr auto native_simd_size_particlereal = 1;

    template<int SIMD_WIDTH = native_simd_size_real>
    using SIMDReal = amrex::Real;

    template<int SIMD_WIDTH = native_simd_size_particlereal>
    using SIMDParticleReal = amrex::ParticleReal;

    // Type that has the same number of int SIMD elements as the SIMDParticleReal type
    template<typename T_ParticleReal = SIMDParticleReal<>>
    using SIMDInt = int;

    // Type that has the same number of IdCpu SIMD elements as the SIMDParticleReal type
    template<typename T_ParticleReal = SIMDParticleReal<>>
    using SIMDIdCpu = std::uint64_t;
#endif

    /// \cond DOXYGEN_IGNORE
    namespace detail {
        struct InternalVectorized {};
    }
    /// \endcond

    /** Mixin Helper Class
     *
     * Use this class as a mixin (base) class to simplify writing functors that support/do not support
     * ParallelForSIMD.
     *
     * Example:
     * ```cpp
     * struct Compute : public Vectorized<>
     * {
     *     template<typename T_Real>
     *     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     *     void operator() (T_Real & AMREX_RESTRICT x) { ... }
     * };
     *
     * // ... call site below
     * {
     *     Compute c;
     *
     *     if constexpr (amrex::simd::is_vectorized<Compute>) {
     *         ParallelForSIMD<...>(np, c).
     *     } else {
     *         ParallelFor(np, c);
     *     }
     * }
     * ```
     */
    template<int SIMD_WIDTH = native_simd_size_real>
    struct
    Vectorized : detail::InternalVectorized
    {
        //! SIMD width in elements
        static constexpr int simd_width = SIMD_WIDTH;
    };

    /** Check if a Functor Class works with amrex::ParallelForSIMD
     *
     * @see amrex::simd::Vectorized
     */
    template<typename T>
    constexpr bool is_vectorized = std::is_base_of_v<detail::InternalVectorized, T>;

    /** Check if a function argument is declared as non-const
     *
     * Use in conjunction with conditional write-back logic from vector registers, e.g.,
     *
     * ```cpp
     * template<typename T_Real>
     * void compute (T_Real & AMREX_RESTRICT x,
     *               T_Real const & AMREX_RESTRICT y) { ... }
     *
     * part_x.copy_from(&m_part_x[i], stdx::element_aligned);
     * part_y.copy_from(&m_part_y[i], stdx::element_aligned);
     *
     * compute(part_x, part_y);
     *
     * if constexpr (is_nth_arg_non_const(compute<double>, 0))
     *     part_x.copy_to(&m_part_x[i], stdx::element_aligned);
     * if constexpr (is_nth_arg_non_const(compute<double>, 1))
     *     part_y.copy_to(&m_part_y[i], stdx::element_aligned);
     * ```
     */
    template<typename R, typename... Args>
    constexpr bool is_nth_arg_non_const (R(*)(Args...), int n)
    {
        constexpr bool val_arr[sizeof...(Args)] {!std::is_const_v<std::remove_reference_t<Args>>...};
        return val_arr[n];
    }
    // same for functors (const/non-const ::operator() members)
    template<typename C, typename R, typename... Args>
    constexpr bool is_nth_arg_non_const (R(C::*)(Args...), int n)
    {
        constexpr bool val_arr[sizeof...(Args)] {!std::is_const_v<std::remove_reference_t<Args>>...};
        return val_arr[n];
    }
    template<typename C, typename R, typename... Args>
    constexpr bool is_nth_arg_non_const (R(C::*)(Args...) const, int n)
    {
        constexpr bool val_arr[sizeof...(Args)] {!std::is_const_v<std::remove_reference_t<Args>>...};
        return val_arr[n];
    }

    /** Load 1D contiguous data from array pointers
     *
     * On GPU and CPU w/o SIMD, this dereferences a 1D array element at the
     * index position i.
     * On CPU with SIMD, this loads a SIMD register at the IndexType::width
     * SIMD-wide index position i.
     *
     * @tparam T data type (amrex::ParticleReal or int or uint64_t)
     * @tparam IndexType int or amrex::SIMDindex<SIMD_WIDTH, int>
     * @param ptr pointer to the array data
     * @param i index or SIMD index
     * @return a reference to the (scalar) data OR a new SIMD variable (register) loaded with data
     */
    template <typename T, typename IndexType>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    decltype(auto) load_1d (T* ptr, IndexType const i)
    {
        if constexpr (std::is_integral_v<IndexType>) {
            return ptr[i];
        } else if constexpr (IndexType::width == 1) {
            return ptr[i.index];
        } else {
#ifdef AMREX_USE_SIMD
            using DataType = stdx::fixed_size_simd<std::decay_t<T>, IndexType::width>;

            // initialize vector register
            // TODO stdx::vector_aligned needs alignment guarantees
            //      https://github.com/AMReX-Codes/amrex/issues/4592
            //      https://en.cppreference.com/w/cpp/experimental/simd/simd/copy_from
            DataType val;
            val.copy_from(&ptr[i.index], stdx::element_aligned);
            return val;

#else
            static_assert(IndexType::width == 1, "SIMD width must be 1 for non-SIMD builds");
            return ptr[i.index];
#endif
        }
    }

    /** Store SIMD register data back to a 1D contiguous array in RAM
     *
     * On GPU and CPU without SIMD, this does nothing because we already
     * modified the (global) RAM directly via pointer.
     *
     * On CPU with SIMD, this performs a conditional writeback of a SIMD register
     * to RAM (index in pointer array), but only if the argument was not passed
     * as const and thus was likely changed.
     *
     * Good optimizing compilers can eliminate writebacks of unchanged values
     * themselves, but we better help a little for robustness. Background:
     * https://github.com/AMReX-Codes/amrex/pull/4520#issuecomment-3064064215
     *
     * @tparam P_Method pointer to the method that used the data (for is_nth_arg_non_const)
     * @tparam N the argument index (for is_nth_arg_non_const)
     * @tparam T data type
     * @tparam IndexType int or SIMD index
     * @tparam ValType the type of the value to store
     * @param val the value to store
     * @param ptr pointer to contiguous 1D array data
     * @param i index or SIMD index
     */
    template <auto P_Method, int N, typename T, typename IndexType, typename ValType>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void store_1d (
        ValType const & AMREX_RESTRICT val,
        T * const AMREX_RESTRICT ptr,
        IndexType const i
    )
    {
        // SIMD uses special vector register types in ValType that need to be copied back to RAM array type T
        if constexpr (!std::is_same_v<ValType, T>) {
            if constexpr (amrex::simd::is_nth_arg_non_const(P_Method, N)) {
#ifdef AMREX_USE_SIMD
                // write back to memory
                // TODO stdx::vector_aligned needs alignment guarantees
                //      https://github.com/AMReX-Codes/amrex/issues/4592
                //      https://en.cppreference.com/w/cpp/experimental/simd/simd/copy_from
                val.copy_to(&ptr[i.index], amrex::simd::stdx::element_aligned);
#else
                amrex::Abort("store_1d: ValType val must alias T ptr data (to make this a no-OP).");
#endif
            }
        }
        amrex::ignore_unused(val, ptr, i);
    }

} // namespace amrex::simd

#endif
