#ifndef AMREX_GPU_LAUNCH_FUNCTS_C_H_
#define AMREX_GPU_LAUNCH_FUNCTS_C_H_
#include <AMReX_Config.H>

namespace amrex {

namespace detail {
    template <typename F, typename N>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, N i)
        noexcept -> decltype(f(0))
    {
        f(i);
    }

    template <typename F, typename N>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, N i)
        noexcept -> decltype(f(0,Gpu::Handler{}))
    {
        f(i,Gpu::Handler{});
    }

    template <typename F>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k)
        noexcept -> decltype(f(0,0,0))
    {
        f(i,j,k);
    }

    template <typename F>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k)
        noexcept -> decltype(f(0,0,0,Gpu::Handler{}))
    {
        f(i,j,k,Gpu::Handler{});
    }

    template <typename F, typename T>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k, T n)
        noexcept -> decltype(f(0,0,0,0))
    {
        f(i,j,k,n);
    }

    template <typename F, typename T>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k, T n)
        noexcept -> decltype(f(0,0,0,0,Gpu::Handler{}))
    {
        f(i,j,k,n,Gpu::Handler{});
    }
}

template<typename T, typename L>
void launch (T const& n, L&& f) noexcept
{
    f(n);
}

template<int MT, typename T, typename L>
void launch (T const& n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    f(n);
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (T n, L&& f) noexcept
{
    for (T i = 0; i < n; ++i) {
        detail::call_f(f,i);
    }
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(n, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    For(n, std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(n, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (T n, L&& f) noexcept
{
    AMREX_PRAGMA_SIMD
    for (T i = 0; i < n; ++i) {
        detail::call_f(f,i);
    }
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(n, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    ParallelFor(n, std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(n, std::forward<L>(f));
}

template <typename L>
void For (Box const& box, L&& f) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        detail::call_f(f,i,j,k);
    }}}
}

template <int MT, typename L>
void For (Box const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box, std::forward<L>(f));
}

template <typename L>
void For (Gpu::KernelInfo const&, Box const& box, L&& f) noexcept
{
    For(box, std::forward<L>(f));
}

template <int MT, typename L>
void For (Gpu::KernelInfo const&, Box const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box, std::forward<L>(f));
}

template <typename L>
void ParallelFor (Box const& box, L&& f) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        detail::call_f(f,i,j,k);
    }}}
}

template <int MT, typename L>
void ParallelFor (Box const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box, std::forward<L>(f));
}

template <typename L>
void ParallelFor (Gpu::KernelInfo const&, Box const& box, L&& f) noexcept
{
    ParallelFor(box, std::forward<L>(f));
}

template <int MT, typename L>
void ParallelFor (Gpu::KernelInfo const&, Box const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Box const& box, T ncomp, L&& f) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (T n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            detail::call_f(f,i,j,k,n);
        }}}
    }
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Box const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box, ncomp, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Gpu::KernelInfo const&, Box const& box, T ncomp, L&& f) noexcept
{
    For(box, ncomp, std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Gpu::KernelInfo const&, Box const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box, ncomp, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (Box const& box, T ncomp, L&& f) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (T n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            detail::call_f(f,i,j,k,n);
        }}}
    }
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (Box const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box, ncomp, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (Gpu::KernelInfo const&, Box const& box, T ncomp, L&& f) noexcept
{
    ParallelFor(box, ncomp, std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (Gpu::KernelInfo const&, Box const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box, ncomp, std::forward<L>(f));
}

template <typename L1, typename L2>
void For (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    For(box1, std::forward<L1>(f1));
    For(box2, std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2>
void For (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1, std::forward<L1>(f1));
    For(box2, std::forward<L2>(f2));
}

template <typename L1, typename L2>
void For (Gpu::KernelInfo const&, Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    For (box1, box2, std::forward<L1>(f1), std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2>
void For (Gpu::KernelInfo const&, Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For (box1, box2, std::forward<L1>(f1), std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void For (Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    For(box1, std::forward<L1>(f1));
    For(box2, std::forward<L2>(f2));
    For(box3, std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3>
void For (Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1, std::forward<L1>(f1));
    For(box2, std::forward<L2>(f2));
    For(box3, std::forward<L3>(f3));
}

template <typename L1, typename L2, typename L3>
void For (Gpu::KernelInfo const&, Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    For(box1, box2, box3, std::forward<L1>(f1), std::forward<L2>(f2), std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3>
void For (Gpu::KernelInfo const&, Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1, box2, box3, std::forward<L1>(f1), std::forward<L2>(f2), std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    For(box1, ncomp1, std::forward<L1>(f1));
    For(box2, ncomp2, std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1, ncomp1, std::forward<L1>(f1));
    For(box2, ncomp2, std::forward<L2>(f2));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void For (Gpu::KernelInfo const&,
          Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void For (Gpu::KernelInfo const&,
          Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2,
          Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    For(box1, ncomp1, std::forward<L1>(f1));
    For(box2, ncomp2, std::forward<L2>(f2));
    For(box3, ncomp3, std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2,
          Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1, ncomp1, std::forward<L1>(f1));
    For(box2, ncomp2, std::forward<L2>(f2));
    For(box3, ncomp3, std::forward<L3>(f3));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void For (Gpu::KernelInfo const&,
          Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2,
          Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void For (Gpu::KernelInfo const&,
          Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2,
          Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <typename L1, typename L2>
void ParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2>
void ParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
}

template <typename L1, typename L2>
void ParallelFor (Gpu::KernelInfo const&, Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(box1,box2,f1,f2);
}

template <int MT, typename L1, typename L2>
void ParallelFor (Gpu::KernelInfo const&, Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,f1,f2);
}

template <typename L1, typename L2, typename L3>
void ParallelFor (Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
    ParallelFor(box3, std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3>
void ParallelFor (Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
    ParallelFor(box3, std::forward<L3>(f3));
}

template <typename L1, typename L2, typename L3>
void ParallelFor (Gpu::KernelInfo const&, Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3>
void ParallelFor (Gpu::KernelInfo const&, Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void ParallelFor (Gpu::KernelInfo const&,
                  Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void ParallelFor (Gpu::KernelInfo const&,
                  Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2,
                  Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
    ParallelFor(box3, ncomp3, std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2,
                  Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
    ParallelFor(box3, ncomp3, std::forward<L3>(f3));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void ParallelFor (Gpu::KernelInfo const&,
                  Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2,
                  Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(box1, ncomp1, std::forward<L1>(f1),
                box2, ncomp2, std::forward<L2>(f2),
                box3, ncomp3, std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void ParallelFor (Gpu::KernelInfo const&,
                  Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2,
                  Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1, ncomp1, std::forward<L1>(f1),
                box2, ncomp2, std::forward<L2>(f2),
                box3, ncomp3, std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (T n, L&& f) noexcept
{
    ParallelFor(n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(n,std::forward<L>(f));
}

template <typename L>
void HostDeviceParallelFor (Box const& box, L&& f) noexcept
{
    ParallelFor(box,std::forward<L>(f));
}

template <int MT, typename L>
void HostDeviceParallelFor (Box const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box,std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (Box const& box, T ncomp, L&& f) noexcept
{
    ParallelFor(box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (Box const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void HostDeviceParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2>
void HostDeviceParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void HostDeviceParallelFor (Box const& box1, Box const& box2, Box const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3>
void HostDeviceParallelFor (Box const& box1, Box const& box2, Box const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2,
                            Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2,
                            Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (T n, L&& f) noexcept
{
    For(n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(n,std::forward<L>(f));
}

template <typename L>
void HostDeviceFor (Box const& box, L&& f) noexcept
{
    For(box,std::forward<L>(f));
}

template <int MT, typename L>
void HostDeviceFor (Box const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box,std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Box const& box, T ncomp, L&& f) noexcept
{
    For(box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Box const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void HostDeviceFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    For(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2>
void HostDeviceFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void HostDeviceFor (Box const& box1, Box const& box2, Box const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    For(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3>
void HostDeviceFor (Box const& box1, Box const& box2, Box const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceFor (Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceFor (Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceFor (Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2,
                    Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceFor (Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2,
                    Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    ParallelFor(n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(n,std::forward<L>(f));
}

template <typename L>
void HostDeviceParallelFor (Gpu::KernelInfo const&, Box const& box, L&& f) noexcept
{
    ParallelFor(box,std::forward<L>(f));
}

template <int MT, typename L>
void HostDeviceParallelFor (Gpu::KernelInfo const&, Box const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box,std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (Gpu::KernelInfo const&, Box const& box, T ncomp, L&& f) noexcept
{
    ParallelFor(box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (Gpu::KernelInfo const&, Box const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void HostDeviceParallelFor (Gpu::KernelInfo const&, Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2>
void HostDeviceParallelFor (Gpu::KernelInfo const&, Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            Box const& box1, Box const& box2, Box const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3>
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            Box const& box1, Box const& box2, Box const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2,
                            Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2,
                            Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    For(n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(n,std::forward<L>(f));
}

template <typename L>
void HostDeviceFor (Gpu::KernelInfo const&, Box const& box, L&& f) noexcept
{
    For(box,std::forward<L>(f));
}

template <int MT, typename L>
void HostDeviceFor (Gpu::KernelInfo const&, Box const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box,std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Gpu::KernelInfo const&, Box const& box, T ncomp, L&& f) noexcept
{
    For(box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Gpu::KernelInfo const&, Box const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void HostDeviceFor (Gpu::KernelInfo const&, Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    For(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2>
void HostDeviceFor (Gpu::KernelInfo const&, Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void HostDeviceFor (Gpu::KernelInfo const&,
                    Box const& box1, Box const& box2, Box const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    For(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3>
void HostDeviceFor (Gpu::KernelInfo const&,
                    Box const& box1, Box const& box2, Box const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceFor (Gpu::KernelInfo const&,
                    Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceFor (Gpu::KernelInfo const&,
                    Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceFor (Gpu::KernelInfo const&,
                    Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2,
                    Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceFor (Gpu::KernelInfo const&,
                    Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2,
                    Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelForRNG (T n, L&& f) noexcept
{
    for (T i = 0; i < n; ++i) {
        f(i,RandomEngine{});
    }
}

template <typename L>
void ParallelForRNG (Box const& box, L&& f) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,RandomEngine{});
    }}}
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelForRNG (Box const& box, T ncomp, L&& f) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (T n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            f(i,j,k,n,RandomEngine{});
        }}}
    }
}

template <typename L>
void single_task (L&& f) noexcept
{
    f();
}

}

#endif
