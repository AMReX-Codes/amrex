#ifndef AMREX_ALG_VEC_UTIL_H_
#define AMREX_ALG_VEC_UTIL_H_
#include <AMReX_Config.H>

#include <AMReX_AlgVector.H>

namespace amrex {

template <class V, class Enable = void> struct IsAlgVector : std::false_type {};
//
template <class V>
struct IsAlgVector<V, std::enable_if_t<std::is_same_v<AlgVector<typename V::value_type,
                                                                typename V::allocator_type>,
                                                      V> > >
    : std::true_type {};

template <typename V1, typename F>
std::enable_if_t<IsAlgVector<std::decay_t<V1> >::value>
ForEach (V1 & x, F const& f)
{
    Long n = x.numLocalRows();
    auto* px = x.data();
    ParallelForOMP(n, [=] AMREX_GPU_DEVICE (Long i) noexcept
    {
        f(px[i]);
    });
}

template <typename V1, typename V2, typename F>
std::enable_if_t<IsAlgVector<std::decay_t<V1> >::value &&
                 IsAlgVector<std::decay_t<V2> >::value>
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

template <typename V1, typename V2, typename V3, typename F>
std::enable_if_t<IsAlgVector<std::decay_t<V1> >::value &&
                 IsAlgVector<std::decay_t<V2> >::value &&
                 IsAlgVector<std::decay_t<V3> >::value>
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

template <typename V1, typename V2, typename V3, typename V4, typename F>
std::enable_if_t<IsAlgVector<std::decay_t<V1> >::value &&
                 IsAlgVector<std::decay_t<V2> >::value &&
                 IsAlgVector<std::decay_t<V3> >::value &&
                 IsAlgVector<std::decay_t<V4> >::value>
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

template <typename V1, typename V2, typename V3, typename V4, typename V5, typename F>
std::enable_if_t<IsAlgVector<std::decay_t<V1> >::value &&
                 IsAlgVector<std::decay_t<V2> >::value &&
                 IsAlgVector<std::decay_t<V3> >::value &&
                 IsAlgVector<std::decay_t<V4> >::value &&
                 IsAlgVector<std::decay_t<V5> >::value>
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

//! Return dot product of two vectors. By default, this returns the global
//! result over all MPI ranks. If `local` is true, it returns the result
//! over the locally stored entries only.
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

//! y = ax + y. For GPU guilds, this function is asynchronous with respect
//! to the host.
template <typename T, typename Allocator>
void Axpy (AlgVector<T,Allocator>& y, T a, AlgVector<T,Allocator> const& x)
{
    ForEach(y, x, [=] AMREX_GPU_DEVICE (T& yi, T const& xi) { yi += a*xi; });
}

//! y = x + a*y. For GPU guilds, this function is asynchronous with respect
//! to the host.
template <typename T, typename Allocator>
void Xpay (AlgVector<T,Allocator>& y, T a, AlgVector<T,Allocator> const& x)
{
    ForEach(y, x, [=] AMREX_GPU_DEVICE (T& yi, T const& xi) { yi = a*yi + xi; });
}

//! y = a*xa + b*xb. For GPU guilds, this function is asynchronous with
//! respect to the host.
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
