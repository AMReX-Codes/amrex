#ifndef AMREX_REDUCE_SIMD_H_
#define AMREX_REDUCE_SIMD_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_OpenMP.H>
#include <AMReX_Reduce.H>  // provides amrex::SIMDindex via AMReX_Gpu.H
#include <AMReX_SIMD.H>

#include <concepts>
#include <cstddef>
#include <limits>
#include <type_traits>

namespace amrex {

/// \cond DOXYGEN_IGNORE
namespace Reduce::detail {

    //! Reduce operations supported by the SIMD reduction path
    template <typename P>
    inline constexpr bool is_simd_reduce_op =
        std::is_same_v<P, ReduceOpSum> ||
        std::is_same_v<P, ReduceOpMin> ||
        std::is_same_v<P, ReduceOpMax>;

#if defined(AMREX_USE_SIMD) && !defined(AMREX_USE_GPU)

    //! Transform a GpuTuple of scalars into a GpuTuple of SIMD registers of width W
    template <int W, typename T>
    struct MakeSIMDTuple;

    template <int W, typename... Ts>
    struct MakeSIMDTuple<W, GpuTuple<Ts...>>
    {
        using type = GpuTuple<amrex::simd::stdx::fixed_size_simd<Ts, W>...>;
    };

    //! Initialize a SIMD accumulator register with the identity element of the reduce operation
    template <typename P, typename T>
    AMREX_FORCE_INLINE
    void simd_init (T& t)
    {
        using E = typename T::value_type;
        if constexpr (std::is_same_v<P, ReduceOpSum>) {
            t = T(E(0));
        } else if constexpr (std::is_same_v<P, ReduceOpMin>) {
            t = T(std::numeric_limits<E>::max());
        } else if constexpr (std::is_same_v<P, ReduceOpMax>) {
            t = T(std::numeric_limits<E>::lowest());
        } else {
            static_assert(is_simd_reduce_op<P>,
                "SIMD reductions support only ReduceOpSum, ReduceOpMin and ReduceOpMax");
        }
    }

    //! Element-wise (per SIMD lane) update of a SIMD accumulator register
    template <typename P, typename T>
    AMREX_FORCE_INLINE
    void simd_local_update (T& d, T const& s)
    {
        namespace stdx = amrex::simd::stdx;
        if constexpr (std::is_same_v<P, ReduceOpSum>) {
            d += s;
        } else if constexpr (std::is_same_v<P, ReduceOpMin>) {
            d = stdx::min(d, s);
        } else if constexpr (std::is_same_v<P, ReduceOpMax>) {
            d = stdx::max(d, s);
        } else {
            static_assert(is_simd_reduce_op<P>,
                "SIMD reductions support only ReduceOpSum, ReduceOpMin and ReduceOpMax");
        }
    }

    //! Horizontal (across SIMD lanes) fold of a SIMD register into a scalar accumulator
    template <typename P, typename T, typename S>
    AMREX_FORCE_INLINE
    void simd_horizontal_update (T& d, S const& s)
    {
        namespace stdx = amrex::simd::stdx;
        if constexpr (std::is_same_v<P, ReduceOpSum>) {
            d += stdx::reduce(s);
        } else if constexpr (std::is_same_v<P, ReduceOpMin>) {
            d = amrex::min(d, stdx::hmin(s));
        } else if constexpr (std::is_same_v<P, ReduceOpMax>) {
            d = amrex::max(d, stdx::hmax(s));
        } else {
            static_assert(is_simd_reduce_op<P>,
                "SIMD reductions support only ReduceOpSum, ReduceOpMin and ReduceOpMax");
        }
    }

    // The recursions below mirror for_each_local / for_each_init in AMReX_Reduce.H.

    template <std::size_t I, typename T, typename P>
    AMREX_FORCE_INLINE
    void for_each_simd_init (T& t)
    {
        simd_init<P>(amrex::get<I>(t));
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_FORCE_INLINE
    void for_each_simd_init (T& t)
    {
        simd_init<P>(amrex::get<I>(t));
        for_each_simd_init<I+1,T,P1,Ps...>(t);
    }

    template <std::size_t I, typename T, typename P>
    AMREX_FORCE_INLINE
    void for_each_simd_local (T& d, T const& s)
    {
        simd_local_update<P>(amrex::get<I>(d), amrex::get<I>(s));
    }

    template <std::size_t I, typename T, typename P, typename P1, typename... Ps>
    AMREX_FORCE_INLINE
    void for_each_simd_local (T& d, T const& s)
    {
        simd_local_update<P>(amrex::get<I>(d), amrex::get<I>(s));
        for_each_simd_local<I+1,T,P1,Ps...>(d, s);
    }

    template <std::size_t I, typename TD, typename TS, typename P>
    AMREX_FORCE_INLINE
    void for_each_horizontal (TD& d, TS const& s)
    {
        simd_horizontal_update<P>(amrex::get<I>(d), amrex::get<I>(s));
    }

    template <std::size_t I, typename TD, typename TS, typename P, typename P1, typename... Ps>
    AMREX_FORCE_INLINE
    void for_each_horizontal (TD& d, TS const& s)
    {
        simd_horizontal_update<P>(amrex::get<I>(d), amrex::get<I>(s));
        for_each_horizontal<I+1,TD,TS,P1,Ps...>(d, s);
    }

#endif // AMREX_USE_SIMD && !AMREX_USE_GPU
}
/// \endcond

#ifndef AMREX_USE_GPU

/**
 * \brief SIMD-vectorized evaluation loop for host (CPU) reductions.
 *
 * This is the CPU/SIMD counterpart of the host \c ReduceOps::eval element loop.
 * It processes \c n elements in chunks of \c WIDTH using per-lane SIMD accumulator
 * registers, handles the remainder with scalar iterations, and folds the result
 * into \c reduce_data exactly like the scalar evaluation loop.
 *
 * The callable \c f is invoked as \c f(amrex::SIMDindex<WIDTH>{i}) in the main loop
 * and as \c f(amrex::SIMDindex<1>{i}) for the remainder. For \c SIMDindex<WIDTH> it
 * must return the SIMD-widened form of the \c ReduceData tuple (each scalar component
 * replaced by a \c fixed_size_simd register of the same element type and width);
 * for \c SIMDindex<1> it must return the scalar \c ReduceData tuple. Kernels written
 * with amrex::simd::load_1d satisfy both forms from a single source.
 *
 * Compiled without SIMD support (\c AMReX_SIMD=OFF), only the scalar loop remains
 * and \c f is only ever called with \c SIMDindex<1>.
 *
 * Only \c ReduceOpSum, \c ReduceOpMin and \c ReduceOpMax are supported.
 *
 * \note SIMD sum reductions reassociate floating-point additions (per-lane partial
 *       sums plus a final horizontal fold), so results can differ from the scalar
 *       evaluation order by rounding.
 *
 * \tparam WIDTH SIMD width in elements
 * \tparam N     index type (integer)
 * \tparam D     the ReduceData type
 * \tparam F     callable, see above
 * \tparam Ps    the reduce operations of the ReduceOps
 *
 * \param n           number of elements to reduce over
 * \param reduce_data where the result is accumulated
 * \param f           callable, see above
 */
template <int WIDTH, std::integral N, typename D, typename F, typename... Ps>
AMREX_ATTRIBUTE_FLATTEN_FOR
void evalReduceSIMD (N n, D& reduce_data, F const& f, ReduceOps<Ps...>& /*reduce_ops*/)
{
    static_assert(WIDTH >= 1, "evalReduceSIMD: WIDTH must be at least 1");
    static_assert((Reduce::detail::is_simd_reduce_op<Ps> && ...),
        "evalReduceSIMD supports only ReduceOpSum, ReduceOpMin and ReduceOpMax");

    using ReduceTuple = typename D::Type;
    ReduceTuple rr;  // scalar accumulator (remainder and horizontal fold target)
    Reduce::detail::for_each_init<0, ReduceTuple, Ps...>(rr);
    N i = 0;
#ifdef AMREX_USE_SIMD
    if constexpr (WIDTH > 1) {
        using SimdTuple = typename Reduce::detail::MakeSIMDTuple<WIDTH, ReduceTuple>::type;
        static_assert(std::is_same_v<SimdTuple,
                          std::decay_t<decltype(f(SIMDindex<WIDTH, N>{0}))>>,
            "evalReduceSIMD: for SIMDindex<WIDTH> the callable must return the "
            "SIMD-widened ReduceData tuple (fixed_size_simd of each component type). "
            "Load particle/array data with amrex::simd::load_1d to get matching types.");
        SimdTuple sr;
        Reduce::detail::for_each_simd_init<0, SimdTuple, Ps...>(sr);
        // vectorized main loop over full SIMD lanes
        for (; i + WIDTH <= n; i += WIDTH) {
            Reduce::detail::for_each_simd_local<0, SimdTuple, Ps...>(sr, f(SIMDindex<WIDTH, N>{i}));
        }
        Reduce::detail::for_each_horizontal<0, ReduceTuple, SimdTuple, Ps...>(rr, sr);
    }
#endif
    // scalar handling of the remainder (or the whole range without SIMD)
    for (; i < n; ++i) {
        Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(rr, f(SIMDindex<1, N>{i}));
    }
    Reduce::detail::for_each_local<0, ReduceTuple, Ps...>(
        reduce_data.reference(OpenMP::get_thread_num()), rr);
}

#endif // !AMREX_USE_GPU
}

#endif
