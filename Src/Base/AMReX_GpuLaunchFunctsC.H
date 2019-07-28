#ifndef AMREX_GPU_LAUNCH_FUNCTS_C_H_
#define AMREX_GPU_LAUNCH_FUNCTS_C_H_

namespace amrex {

template<typename T, typename L>
AMREX_FORCE_INLINE
void launch (T const& n, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    f(n);
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
void For (T n, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    for (T i = 0; i < n; ++i) {
        f(i);
    }
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
void ParallelFor (T n, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    AMREX_PRAGMA_SIMD
    for (T i = 0; i < n; ++i) {
        f(i);
    }
}

template <typename L>
AMREX_FORCE_INLINE
void For (Box const& box, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <typename L>
AMREX_FORCE_INLINE
void ParallelFor (Box const& box, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
void For (Box const& box, T ncomp, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (T n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            f(i,j,k,n);
        }}}
    }
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
AMREX_FORCE_INLINE
void ParallelFor (Box const& box, T ncomp, L&& f, std::size_t shared_mem_bytes=0) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (T n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            f(i,j,k,n);
        }}}
    }
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void FabReduce (Box const& box, N ncomp, T const& init_val,
                L1&& f1, L2&& f2, std::size_t shared_mem_bytes=0) noexcept
{
    auto r = init_val;
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (N n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            f1(i,j,k,n,&r);
        }}}
    }
    f2(r);
}

template <typename T, typename L1, typename L2>
void FabReduce (Box const& box, T const& init_val,
                L1&& f1, L2&& f2, std::size_t shared_mem_bytes=0) noexcept
{
    auto r = init_val;
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f1(i,j,k,&r);
    }}}
    f2(r);
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void VecReduce (N n, T const& init_val,
                L1&& f1, L2&& f2, std::size_t shared_mem_bytes=0) noexcept
{
    auto r = init_val;
    AMREX_PRAGMA_SIMD
    for (N i = 0; i < n; ++i) {
        f1(i,&r);
    }
    f2(r);
}

namespace Reduce {

template <typename N, typename T, typename U, typename BOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Sum (N n, U const* v, T init_val, BOP bop)
{
    T sum = init_val;
    T* dp = &sum;
    amrex::VecReduce(n, init_val,
    [=] (N i, T* r) noexcept
    {
        *r = bop(*r, v[i]);
    },
    [=] (T r) noexcept
    {
        Gpu::ReduceSum(dp, r);
    });
    return sum;
}

template <typename N, typename T, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Sum (N n, T const* v, T init_val = 0)
{
    return Reduce::Sum(n, v, init_val, amrex::Plus<T>());
}

template <typename N, typename T, typename U, typename BOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Min (N n, U const* v, T init_val, BOP bop)
{
    T mn = init_val;
    T* dp = &mn;
    amrex::VecReduce(n, init_val,
    [=] AMREX_GPU_DEVICE (N i, T* r) noexcept
    {
        *r = bop(*r, v[i]);
    },
    [=] AMREX_GPU_DEVICE (T r) noexcept
    {
        Gpu::ReduceMin(dp, r);
    });
    return mn;
}

template <typename N, typename T, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Min (N n, T const* v, T init_val = std::numeric_limits<T>::max())
{
    return Reduce::Min(n, v, init_val, amrex::Less<T>());
}

template <typename N, typename T, typename U, typename BOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Max (N n, U const* v, T init_val, BOP bop)
{
    T mx = init_val;
    T* dp = &mx;
    amrex::VecReduce(n, init_val,
    [=] AMREX_GPU_DEVICE (N i, T* r) noexcept
    {
        *r = bop(*r, v[i]);
    },
    [=] AMREX_GPU_DEVICE (T r) noexcept
    {
        Gpu::ReduceMax(dp, r);
    });
    return mx;
}

template <typename N, typename T, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
T Max (N n, T const* v, T init_val = std::numeric_limits<T>::lowest())
{
    return Reduce::Max(n, v, init_val, amrex::Greater<T>());
}

template <typename T, typename N, typename U, typename MINOP, typename MAXOP, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
std::pair<T,T> MinMax (N n, U const* v, MINOP minop, MAXOP maxop)
{
    Array<T,2> hv{std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest()};
    T* dp = hv.data();
    typedef GpuArray<T,2> Real2;
    amrex::VecReduce(n, Real2{hv[0],hv[1]},
    [=] AMREX_GPU_DEVICE (N i, Real2* r) noexcept
    {
        (*r)[0] = minop((*r)[0], v[i]);
        (*r)[1] = maxop((*r)[1], v[i]);
    },
    [=] AMREX_GPU_DEVICE (Real2 const& r) noexcept
    {
        Gpu::ReduceMin(dp  , r[0]);
        Gpu::ReduceMax(dp+1, r[1]);
    });
    return std::make_pair(hv[0],hv[1]);
}

template <typename T, typename N, typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
std::pair<T,T> MinMax (N n, T const* v)
{
    return Reduce::MinMax<T>(n, v, amrex::Less<T>(), amrex::Greater<T>());
}

}

}

#endif
