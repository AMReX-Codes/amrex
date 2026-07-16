#ifndef AMREX_CONSTEXPR_FOR_H_
#define AMREX_CONSTEXPR_FOR_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_ConstexprFor.H
 * \brief Compile-time unrolled loop utility.
 */

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>

#include <type_traits>

namespace amrex {

/**
 * \ingroup amrex_utilities
 * \brief Compile-time unrolled loop from \p I (inclusive) to \p N (exclusive).
 *
 * Approximates the compile-time unrolling of:
 * \code{.cpp}
 *     for (auto i = I; i < N; ++i) { f(i); }
 * \endcode
 *
 * The body \p f is a callable (typically a lambda) that accepts a single
 * integral argument. This can be used inside GPU kernels.
 *
 * \note Based on https://artificial-mind.net/blog/2020/10/31/constexpr-for
 *
 * \tparam I Starting index (compile-time constant).
 * \tparam N One-past-end index (compile-time constant).
 * \tparam F Callable type; must accept an integral argument.
 * \param  f Callable invoked for each index from \p I to \p N-1.
 */
template<auto I, auto N, class F>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr void constexpr_for (F const& f)
{
    if constexpr (I < N) {
        f(std::integral_constant<decltype(I), I>());
        constexpr_for<I+1, N>(f);
    }
}

}

#endif
