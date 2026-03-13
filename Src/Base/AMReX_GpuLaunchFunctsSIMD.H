#ifndef AMREX_GPU_LAUNCH_FUNCTS_SIMD_H_
#define AMREX_GPU_LAUNCH_FUNCTS_SIMD_H_

#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_SIMD.H>

#include <type_traits>


namespace amrex
{

/** Helper type to store/access the SIMD width in ParallelForSIMD lambdas
 *
 * Use instead of int as the running index i. Used to pass the
 * SIMD WIDTH as compile-time meta-data into a called function/method.
 *
 * @tparam WIDTH SIMD width in elements
 * @tparam N index type (integer)
 */
template<int WIDTH, class N=int>
struct SIMDindex
{
    /** SIMD width in elements */
    static constexpr int width = WIDTH;

    /** The linear loop index of ParallelFor(SIMD) */
    N index = 0;
};

/** ParallelFor with a SIMD Width (in elements)
 *
 * SIMD load/Write-back operations need to be performed before/after calling this.
 *
 * @tparam WIDTH SIMD width in elements
 * @tparam N index type (integer)
 * @tparam L function/functor to call per SIMD set of elements
 */
template <int WIDTH, typename N, typename L, typename M=std::enable_if_t<std::is_integral_v<N>> >
AMREX_ATTRIBUTE_FLATTEN_FOR
void ParallelForSIMD (N n, L const& f) noexcept
{
    N i = 0;
    // vectorize full lanes
    for (; i + WIDTH <= n; i+=WIDTH) {
        f(SIMDindex<WIDTH, N>{i});
    }
    // scalar handling of the remainder
    // note: we could make the remainder calls faster, by repeatedly
    //       decreasing the SIMD width by 2 until we reach 1
    for (; i < n; ++i) {
        f(SIMDindex<1, N>{i});
    }
}

/** ParallelFor with a SIMD Width (derived from a type w/ fallback to non-SIMD)
 *
 * SIMD load/Write-back operations need to be performed before/after calling this.
 *
 * @tparam T a type that implements amrex::simd::Vectorized (or not)
 * @tparam N index type (integer)
 * @tparam L function/functor to call per SIMD set of elements
 */
template <typename T, typename N, typename L, typename M=std::enable_if_t<std::is_integral_v<N>> >
AMREX_ATTRIBUTE_FLATTEN_FOR
void ParallelForSIMD (N n, L && f) noexcept
{
#ifdef AMREX_USE_SIMD
    if constexpr (amrex::simd::is_vectorized<T>) {
        amrex::ParallelForSIMD<T::simd_width>(n, std::forward<L>(f));
    } else
#endif
    {
        amrex::ParallelFor(n, std::forward<L>(f));
    }
}

} // namespace amrex

#endif
