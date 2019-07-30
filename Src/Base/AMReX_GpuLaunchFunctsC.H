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

}

#endif
