#ifndef AMREX_GPU_LAUNCH_FUNCTS_C_H_
#define AMREX_GPU_LAUNCH_FUNCTS_C_H_

namespace amrex {

template<typename T, typename L>
void launch (T const& n, L&& f, std::size_t /*shared_mem_bytes*/=0) noexcept
{
    f(n);
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void For (T n, L&& f, std::size_t /*shared_mem_bytes*/=0) noexcept
{
    for (T i = 0; i < n; ++i) {
        f(i);
    }
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void ParallelFor (T n, L&& f, std::size_t /*shared_mem_bytes*/=0) noexcept
{
    AMREX_PRAGMA_SIMD
    for (T i = 0; i < n; ++i) {
        f(i);
    }
}

template <typename L>
void For (Box const& box, L&& f, std::size_t /*shared_mem_bytes*/=0) noexcept
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
void ParallelFor (Box const& box, L&& f, std::size_t /*shared_mem_bytes*/=0) noexcept
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
void For (Box const& box, T ncomp, L&& f, std::size_t /*shared_mem_bytes*/=0) noexcept
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
void ParallelFor (Box const& box, T ncomp, L&& f, std::size_t /*shared_mem_bytes*/=0) noexcept
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

template <typename L1, typename L2>
void For (Box const& box1, Box const& box2, L1&& f1, L2&& f2,
          std::size_t /*shared_mem_bytes*/=0) noexcept
{
    For(box1, std::forward<L1>(f1));
    For(box2, std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void For (Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3,
          std::size_t /*shared_mem_bytes*/=0) noexcept
{
    For(box1, std::forward<L1>(f1));
    For(box2, std::forward<L2>(f2));
    For(box3, std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2,
          std::size_t /*shared_mem_bytes*/=0) noexcept
{
    For(box1, ncomp1, std::forward<L1>(f1));
    For(box2, ncomp2, std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value>,
          typename M3=amrex::EnableIf_t<std::is_integral<T3>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2,
          Box const& box3, T3 ncomp3, L3&& f3,
          std::size_t /*shared_mem_bytes*/=0) noexcept
{
    For(box1, ncomp1, std::forward<L1>(f1));
    For(box2, ncomp2, std::forward<L2>(f2));
    For(box3, ncomp3, std::forward<L3>(f3));
}

template <typename L1, typename L2>
void ParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2,
                  std::size_t /*shared_mem_bytes*/=0) noexcept
{
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void ParallelFor (Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3,
                  std::size_t /*shared_mem_bytes*/=0) noexcept
{
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
    ParallelFor(box3, std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2,
                  std::size_t /*shared_mem_bytes*/=0) noexcept
{
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value>,
          typename M3=amrex::EnableIf_t<std::is_integral<T3>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2,
                  Box const& box3, T3 ncomp3, L3&& f3,
                  std::size_t /*shared_mem_bytes*/=0) noexcept
{
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
    ParallelFor(box3, ncomp3, std::forward<L3>(f3));
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void FabReduce (Box const& box, N ncomp, T const& init_val,
                L1&& f1, L2&& f2, std::size_t /*shared_mem_bytes*/=0) noexcept
{
    auto r = init_val;
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (N n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            f1(i,j,k,n,&r);
        }}}
    }
    f2(r);
}

template <typename T, typename L1, typename L2>
void FabReduce (Box const& box, T const& init_val,
                L1&& f1, L2&& f2, std::size_t /*shared_mem_bytes*/=0) noexcept
{
    auto r = init_val;
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f1(i,j,k,&r);
    }}}
    f2(r);
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void VecReduce (N n, T const& init_val,
                L1&& f1, L2&& f2, std::size_t /*shared_mem_bytes*/=0) noexcept
{
    auto r = init_val;
    for (N i = 0; i < n; ++i) {
        f1(i,&r);
    }
    f2(r);
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (T n, L&& f) noexcept
{
    ParallelFor(n,std::forward<L>(f));
}

template <typename L>
void HostDeviceParallelFor (Box const& box, L&& f) noexcept
{
    ParallelFor(box,std::forward<L>(f));
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (Box const& box, T ncomp, L&& f) noexcept
{
    ParallelFor(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void HostDeviceParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void HostDeviceParallelFor (Box const& box1, Box const& box2, Box const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value> >
void HostDeviceParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value>,
          typename M3=amrex::EnableIf_t<std::is_integral<T3>::value> >
void HostDeviceParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2,
                            Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void HostDeviceFor (T n, L&& f) noexcept
{
    For(n,std::forward<L>(f));
}

template <typename L>
void HostDeviceFor (Box const& box, L&& f) noexcept
{
    For(box,std::forward<L>(f));
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void HostDeviceFor (Box const& box, T ncomp, L&& f) noexcept
{
    For(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void HostDeviceFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    For(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void HostDeviceFor (Box const& box1, Box const& box2, Box const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    For(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value> >
void HostDeviceFor (Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value>,
          typename M3=amrex::EnableIf_t<std::is_integral<T3>::value> >
void HostDeviceFor (Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2,
                    Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

}

#endif
