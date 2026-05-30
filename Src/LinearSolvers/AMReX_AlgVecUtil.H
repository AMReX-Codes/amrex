#ifndef AMREX_ALG_VEC_UTIL_H_
#define AMREX_ALG_VEC_UTIL_H_
#include <AMReX_Config.H>

#include <AMReX_AlgVector.H>

namespace amrex {

/**
 * \file AMReX_AlgVecUtil.H
 *
 * Device-friendly helpers for element-wise operations on AlgVector
 * specializations (ForEach, BLAS-like routines, etc.).
 */

template <class V> struct IsAlgVector : std::false_type {};
//
template <class V>
requires (std::same_as<AlgVector<typename V::value_type,
                                 typename V::allocator_type>,
                       V>)
struct IsAlgVector<V> : std::true_type {};

/**
 * \brief Apply \p f to every element of \p x.
 *
 * \param x AlgVector whose entries are visited in place.
 * \param f Callable with signature `void(Value&)` that will be executed on the
 *          device/host for each local row.
 */
template <typename V1, typename F>
requires (IsAlgVector<std::decay_t<V1> >::value)
void
ForEach (V1 & x, F const& f)
{
    Long n = x.numLocalRows();
    auto* px = x.data();
    ParallelForOMP(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        f(px[i]);
    });
}

/**
 * \brief Apply \p f to paired entries from \p x and \p y.
 *
 * \param x First AlgVector; entries may be modified by \p f.
 * \param y Second AlgVector; entries may be modified by \p f.
 * \param f Callable with signature `void(X&, Y&)`.
 */
template <typename V1, typename V2, typename F>
requires (IsAlgVector<std::decay_t<V1> >::value && IsAlgVector<std::decay_t<V2> >::value)
void
ForEach (V1 & x, V2 & y, F const& f)
{
    AMREX_ASSERT(x.numLocalRows() == y.numLocalRows());
    Long n = x.numLocalRows();
    auto* AMREX_RESTRICT px = x.data();
    auto* AMREX_RESTRICT py = y.data();
    ParallelForOMP(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        f(px[i], py[i]);
    });
}

/**
 * \brief Apply \p f to triplets drawn from \p x, \p y, and \p z.
 *
 * \param x First AlgVector; entries may be modified by \p f.
 * \param y Second AlgVector; entries may be modified by \p f.
 * \param z Third AlgVector; entries may be modified by \p f.
 * \param f Callable with signature `void(X&, Y&, Z&)`.
 */
template <typename V1, typename V2, typename V3, typename F>
requires (IsAlgVector<std::decay_t<V1> >::value && IsAlgVector<std::decay_t<V2> >::value && IsAlgVector<std::decay_t<V3> >::value)
void
ForEach (V1 & x, V2 & y, V3 & z, F const& f)
{
    AMREX_ASSERT(x.numLocalRows() == y.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == z.numLocalRows());
    Long n = x.numLocalRows();
    auto* AMREX_RESTRICT px = x.data();
    auto* AMREX_RESTRICT py = y.data();
    auto* AMREX_RESTRICT pz = z.data();
    ParallelForOMP(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        f(px[i], py[i], pz[i]);
    });
}

/**
 * \brief Apply \p f to tuples drawn from four AlgVectors.
 *
 * \param x First AlgVector; entries may be modified by \p f.
 * \param y Second AlgVector; entries may be modified by \p f.
 * \param z Third AlgVector; entries may be modified by \p f.
 * \param a Fourth AlgVector; entries may be modified by \p f.
 * \param f Callable with signature `void(X&, Y&, Z&, A&)`.
 */
template <typename V1, typename V2, typename V3, typename V4, typename F>
requires (IsAlgVector<std::decay_t<V1> >::value && IsAlgVector<std::decay_t<V2> >::value && IsAlgVector<std::decay_t<V3> >::value && IsAlgVector<std::decay_t<V4> >::value)
void
ForEach (V1 & x, V2 & y, V3 & z, V4 & a, F const& f)
{
    AMREX_ASSERT(x.numLocalRows() == y.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == z.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == a.numLocalRows());
    Long n = x.numLocalRows();
    auto* AMREX_RESTRICT px = x.data();
    auto* AMREX_RESTRICT py = y.data();
    auto* AMREX_RESTRICT pz = z.data();
    auto* AMREX_RESTRICT pa = a.data();
    ParallelForOMP(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        f(px[i], py[i], pz[i], pa[i]);
    });
}

/**
 * \brief Apply \p f to tuples spanning five AlgVectors.
 *
 * \param x First AlgVector; entries may be modified by \p f.
 * \param y Second AlgVector; entries may be modified by \p f.
 * \param z Third AlgVector; entries may be modified by \p f.
 * \param a Fourth AlgVector; entries may be modified by \p f.
 * \param b Fifth AlgVector; entries may be modified by \p f.
 * \param f Callable with signature `void(X&, Y&, Z&, A&, B&)`.
 */
template <typename V1, typename V2, typename V3, typename V4, typename V5, typename F>
requires (IsAlgVector<std::decay_t<V1> >::value &&
          IsAlgVector<std::decay_t<V2> >::value &&
          IsAlgVector<std::decay_t<V3> >::value &&
          IsAlgVector<std::decay_t<V4> >::value &&
          IsAlgVector<std::decay_t<V5> >::value)
void
ForEach (V1 & x, V2 & y, V3 & z, V4 & a, V5 & b, F const& f)
{
    AMREX_ASSERT(x.numLocalRows() == y.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == z.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == a.numLocalRows());
    AMREX_ASSERT(x.numLocalRows() == b.numLocalRows());
    Long n = x.numLocalRows();
    auto* AMREX_RESTRICT px = x.data();
    auto* AMREX_RESTRICT py = y.data();
    auto* AMREX_RESTRICT pz = z.data();
    auto* AMREX_RESTRICT pa = a.data();
    auto* AMREX_RESTRICT pb = b.data();
    ParallelForOMP(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        f(px[i], py[i], pz[i], pa[i], pb[i]);
    });
}

/**
 * \brief Dot product of two AlgVectors with optional local-only reduction.
 *
 * \param x First vector.
 * \param y Second vector.
 * \param local Set true to skip the MPI reduction and return only this rank's
 *              contribution.
 * \return Dot product of \p x and \p y (global unless \p local is true).
 */
template <typename T, typename Allocator>
T Dot (AlgVector<T,Allocator> const& x, AlgVector<T,Allocator> const& y, bool local = false)
{
    AMREX_ASSERT(x.numLocalRows() == y.numLocalRows());
    Long n = x.numLocalRows();
    auto const* px = x.data();
    auto const* py = y.data();
    T r = Reduce::Sum<T>(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
                         {
                             return px[i] * py[i];
                         });
    if (!local) {
        ParallelAllReduce::Sum(r, ParallelContext::CommunicatorSub());
    }
    return r;
}

//! y = ax + y. For GPU builds this is asynchronous with respect to the host.
template <typename T, typename Allocator>
void Axpy (AlgVector<T,Allocator>& y, T a, AlgVector<T,Allocator> const& x)
{
    ForEach(y, x, [=] AMREX_GPU_DEVICE (T& yi, T const& xi) { yi += a*xi; });
}

//! y = x + a*y. For GPU builds this is asynchronous with respect to the host.
template <typename T, typename Allocator>
void Xpay (AlgVector<T,Allocator>& y, T a, AlgVector<T,Allocator> const& x)
{
    ForEach(y, x, [=] AMREX_GPU_DEVICE (T& yi, T const& xi) { yi = a*yi + xi; });
}

//! y = a*xa + b*xb. For GPU builds this is asynchronous with respect to the host.
template <typename T, typename Allocator>
void LinComb (AlgVector<T,Allocator>& y, T a, AlgVector<T,Allocator> const& xa,
              T b, AlgVector<T,Allocator> const& xb)
{
    ForEach(y, xa, xb, [=] AMREX_GPU_DEVICE (T& yi, T const& xai, T const& xbi) {
                           yi = a*xai + b*xbi;
                       });
}

}

#endif
