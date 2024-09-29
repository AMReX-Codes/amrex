#ifndef AMREX_GPU_LAUNCH_FUNCTS_C_H_
#define AMREX_GPU_LAUNCH_FUNCTS_C_H_
#include <AMReX_Config.H>

namespace amrex {

namespace detail {

    // call_f_scalar_handler

    template <typename F, typename N>
    AMREX_FORCE_INLINE
    auto call_f_scalar_handler (F const& f, N i)
        noexcept -> decltype(f(0))
    {
        f(i);
    }

    template <typename F, typename N>
    AMREX_FORCE_INLINE
    auto call_f_scalar_handler (F const& f, N i)
        noexcept -> decltype(f(0,Gpu::Handler{}))
    {
        f(i, Gpu::Handler{});
    }

    // call_f_intvect_inner

    template <typename F, std::size_t...Ns, class...Args>
    AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<1> iv, Args...args)
        noexcept -> decltype(f(0, 0, 0, args...))
    {
        f(iv[0], 0, 0, args...);
    }

    template <typename F, std::size_t...Ns, class...Args>
    AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<2> iv, Args...args)
        noexcept -> decltype(f(0, 0, 0, args...))
    {
        f(iv[0], iv[1], 0, args...);
    }

    template <typename F, int dim, std::size_t...Ns, class...Args>
    AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<dim> iv, Args...args)
        noexcept -> decltype(f(iv, args...))
    {
        f(iv, args...);
    }

    template <typename F, int dim, std::size_t...Ns, class...Args>
    AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<dim> iv, Args...args)
        noexcept -> decltype(f(iv[Ns]..., args...))
    {
        f(iv[Ns]..., args...);
    }

    // call_f_intvect_engine

    template <typename F, int dim>
    AMREX_FORCE_INLINE
    auto call_f_intvect_engine (F const& f, IntVectND<dim> iv, RandomEngine engine)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, engine))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, engine);
    }

    // call_f_intvect_handler

    template <typename F, int dim>
    AMREX_FORCE_INLINE
    auto call_f_intvect_handler (F const& f, IntVectND<dim> iv)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv);
    }

    template <typename F, int dim>
    AMREX_FORCE_INLINE
    auto call_f_intvect_handler (F const& f, IntVectND<dim> iv)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, Gpu::Handler{}))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, Gpu::Handler{});
    }

    // call_f_intvect_ncomp_engine

    template <typename F, typename T, int dim>
    AMREX_FORCE_INLINE
    auto call_f_intvect_ncomp_engine (F const& f, IntVectND<dim> iv, T n, RandomEngine engine)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n, engine))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n, engine);
    }

    // call_f_intvect_ncomp_handler

    template <typename F, typename T, int dim>
    AMREX_FORCE_INLINE
    auto call_f_intvect_ncomp_handler (F const& f, IntVectND<dim> iv, T n)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n);
    }

    template <typename F, typename T, int dim>
    AMREX_FORCE_INLINE
    auto call_f_intvect_ncomp_handler (F const& f, IntVectND<dim> iv, T n)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n, Gpu::Handler{}))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n, Gpu::Handler{});
    }

}

template<typename T, typename L>
void launch (T const& n, L&& f) noexcept
{
    std::forward<L>(f)(n);
}

template<int MT, typename T, typename L>
void launch (T const& n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    std::forward<L>(f)(n);
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
AMREX_ATTRIBUTE_FLATTEN_FOR
void For (T n, L const& f) noexcept
{
    for (T i = 0; i < n; ++i) {
        detail::call_f_scalar_handler(f,i);
    }
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void For (T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(n, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void For (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    For(n, std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void For (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(n, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
AMREX_ATTRIBUTE_FLATTEN_FOR
void ParallelFor (T n, L const& f) noexcept
{
    AMREX_PRAGMA_SIMD
    for (T i = 0; i < n; ++i) {
        detail::call_f_scalar_handler(f,i);
    }
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void ParallelFor (T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(n, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void ParallelFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    ParallelFor(n, std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void ParallelFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(n, std::forward<L>(f));
}

namespace detail {

template <int idim, typename L, int dim>
AMREX_FORCE_INLINE
void For_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv) noexcept
{
    if constexpr (idim == 1) {
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_handler(f,iv);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_handler(f,iv);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_handler(f,iv);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            For_impND<idim-1>(f, lo, hi, iv);
        }
    }
}

}

template <typename L, int dim>
AMREX_ATTRIBUTE_FLATTEN_FOR
void For (BoxND<dim> const& box, L const& f) noexcept
{
    const auto lo = amrex::lbound_iv(box);
    const auto hi = amrex::ubound_iv(box);
    IntVectND<dim> iv;
    detail::For_impND<dim>(f, lo, hi, iv);
}

template <int MT, typename L, int dim>
void For (BoxND<dim> const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box, std::forward<L>(f));
}

template <typename L, int dim>
void For (Gpu::KernelInfo const&, BoxND<dim> const& box, L&& f) noexcept
{
    For(box, std::forward<L>(f));
}

template <int MT, typename L, int dim>
void For (Gpu::KernelInfo const&, BoxND<dim> const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box, std::forward<L>(f));
}

namespace detail {

template <int idim, typename L, int dim>
AMREX_FORCE_INLINE
void ParallelFor_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv) noexcept
{
    if constexpr (idim == 1) {
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_handler(f,iv);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_handler(f,iv);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_handler(f,iv);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            ParallelFor_impND<idim-1>(f, lo, hi, iv);
        }
    }
}

}

template <typename L, int dim>
AMREX_ATTRIBUTE_FLATTEN_FOR
void ParallelFor (BoxND<dim> const& box, L const& f) noexcept
{
    const auto lo = amrex::lbound_iv(box);
    const auto hi = amrex::ubound_iv(box);
    IntVectND<dim> iv;
    detail::ParallelFor_impND<dim>(f, lo, hi, iv);
}

template <int MT, typename L, int dim>
void ParallelFor (BoxND<dim> const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box, std::forward<L>(f));
}

template <typename L, int dim>
void ParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box, L&& f) noexcept
{
    ParallelFor(box, std::forward<L>(f));
}

template <int MT, typename L, int dim>
void ParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box, std::forward<L>(f));
}

namespace detail {

template <int idim, typename L, typename T, int dim>
AMREX_FORCE_INLINE
void For_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv, T n) noexcept
{
    if constexpr (idim == 1) {
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_handler(f,iv,n);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_handler(f,iv,n);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_handler(f,iv,n);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            For_impND<idim-1>(f, lo, hi, iv, n);
        }
    }
}

}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
AMREX_ATTRIBUTE_FLATTEN_FOR
void For (BoxND<dim> const& box, T ncomp, L const& f) noexcept
{
    const auto lo = amrex::lbound_iv(box);
    const auto hi = amrex::ubound_iv(box);
    IntVectND<dim> iv;
    for (T n = 0; n < ncomp; ++n) {
        detail::For_impND<dim>(f, lo, hi, iv, n);
    }
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void For (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box, ncomp, std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void For (Gpu::KernelInfo const&, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    For(box, ncomp, std::forward<L>(f));
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void For (Gpu::KernelInfo const&, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box, ncomp, std::forward<L>(f));
}

namespace detail {

template <int idim, typename L, typename T, int dim>
AMREX_FORCE_INLINE
void ParallelFor_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv, T n) noexcept
{
    if constexpr (idim == 1) {
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_handler(f,iv,n);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_handler(f,iv,n);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_handler(f,iv,n);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            ParallelFor_impND<idim-1>(f, lo, hi, iv, n);
        }
    }
}

}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
AMREX_ATTRIBUTE_FLATTEN_FOR
void ParallelFor (BoxND<dim> const& box, T ncomp, L const& f) noexcept
{
    const auto lo = amrex::lbound_iv(box);
    const auto hi = amrex::ubound_iv(box);
    IntVectND<dim> iv;
    for (T n = 0; n < ncomp; ++n) {
        detail::ParallelFor_impND<dim>(f, lo, hi, iv, n);
    }
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void ParallelFor (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box, ncomp, std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void ParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    ParallelFor(box, ncomp, std::forward<L>(f));
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void ParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box, ncomp, std::forward<L>(f));
}

template <typename L1, typename L2, int dim>
void For (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    For(box1, std::forward<L1>(f1));
    For(box2, std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void For (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1, std::forward<L1>(f1));
    For(box2, std::forward<L2>(f2));
}

template <typename L1, typename L2, int dim>
void For (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    For (box1, box2, std::forward<L1>(f1), std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void For (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For (box1, box2, std::forward<L1>(f1), std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
void For (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    For(box1, std::forward<L1>(f1));
    For(box2, std::forward<L2>(f2));
    For(box3, std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void For (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1, std::forward<L1>(f1));
    For(box2, std::forward<L2>(f2));
    For(box3, std::forward<L3>(f3));
}

template <typename L1, typename L2, typename L3, int dim>
void For (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    For(box1, box2, box3, std::forward<L1>(f1), std::forward<L2>(f2), std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void For (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1, box2, box3, std::forward<L1>(f1), std::forward<L2>(f2), std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void For (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    For(box1, ncomp1, std::forward<L1>(f1));
    For(box2, ncomp2, std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void For (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1, ncomp1, std::forward<L1>(f1));
    For(box2, ncomp2, std::forward<L2>(f2));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void For (Gpu::KernelInfo const&,
          BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void For (Gpu::KernelInfo const&,
          BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void For (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
          BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    For(box1, ncomp1, std::forward<L1>(f1));
    For(box2, ncomp2, std::forward<L2>(f2));
    For(box3, ncomp3, std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void For (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
          BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1, ncomp1, std::forward<L1>(f1));
    For(box2, ncomp2, std::forward<L2>(f2));
    For(box3, ncomp3, std::forward<L3>(f3));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void For (Gpu::KernelInfo const&,
          BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
          BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void For (Gpu::KernelInfo const&,
          BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
          BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <typename L1, typename L2, int dim>
void ParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void ParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
}

template <typename L1, typename L2, int dim>
void ParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void ParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
void ParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
    ParallelFor(box3, std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void ParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
    ParallelFor(box3, std::forward<L3>(f3));
}

template <typename L1, typename L2, typename L3, int dim>
void ParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void ParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void ParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void ParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void ParallelFor (Gpu::KernelInfo const&,
                  BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void ParallelFor (Gpu::KernelInfo const&,
                  BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void ParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                  BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
    ParallelFor(box3, ncomp3, std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void ParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                  BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
    ParallelFor(box3, ncomp3, std::forward<L3>(f3));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void ParallelFor (Gpu::KernelInfo const&,
                  BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                  BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(box1, ncomp1, std::forward<L1>(f1),
                box2, ncomp2, std::forward<L2>(f2),
                box3, ncomp3, std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void ParallelFor (Gpu::KernelInfo const&,
                  BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                  BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1, ncomp1, std::forward<L1>(f1),
                box2, ncomp2, std::forward<L2>(f2),
                box3, ncomp3, std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceParallelFor (T n, L&& f) noexcept
{
    ParallelFor(n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceParallelFor (T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(n,std::forward<L>(f));
}

template <typename L, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box, L&& f) noexcept
{
    ParallelFor(box,std::forward<L>(f));
}

template <int MT, typename L, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box,std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceParallelFor (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    ParallelFor(box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceParallelFor (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void HostDeviceParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void HostDeviceParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void HostDeviceParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                            BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void HostDeviceParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                            BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceFor (T n, L&& f) noexcept
{
    For(n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceFor (T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(n,std::forward<L>(f));
}

template <typename L, int dim>
void HostDeviceFor (BoxND<dim> const& box, L&& f) noexcept
{
    For(box,std::forward<L>(f));
}

template <int MT, typename L, int dim>
void HostDeviceFor (BoxND<dim> const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box,std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceFor (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    For(box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, int dim, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceFor (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2, int dim>
void HostDeviceFor (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    For(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void HostDeviceFor (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
void HostDeviceFor (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    For(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void HostDeviceFor (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void HostDeviceFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void HostDeviceFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void HostDeviceFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                    BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void HostDeviceFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                    BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceParallelFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    ParallelFor(n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceParallelFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(n,std::forward<L>(f));
}

template <typename L, int dim>
void HostDeviceParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box, L&& f) noexcept
{
    ParallelFor(box,std::forward<L>(f));
}

template <int MT, typename L, int dim>
void HostDeviceParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box,std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    ParallelFor(box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2, int dim>
void HostDeviceParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void HostDeviceParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                            BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void HostDeviceParallelFor (Gpu::KernelInfo const&,
                            BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                            BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    For(n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceFor (Gpu::KernelInfo const&, T n, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(n,std::forward<L>(f));
}

template <typename L, int dim>
void HostDeviceFor (Gpu::KernelInfo const&, BoxND<dim> const& box, L&& f) noexcept
{
    For(box,std::forward<L>(f));
}

template <int MT, typename L, int dim>
void HostDeviceFor (Gpu::KernelInfo const&, BoxND<dim> const& box, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box,std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceFor (Gpu::KernelInfo const&, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    For(box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
void HostDeviceFor (Gpu::KernelInfo const&, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    amrex::ignore_unused(MT);
    For(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2, int dim>
void HostDeviceFor (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    For(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void HostDeviceFor (Gpu::KernelInfo const&, BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
void HostDeviceFor (Gpu::KernelInfo const&,
                    BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    For(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void HostDeviceFor (Gpu::KernelInfo const&,
                    BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void HostDeviceFor (Gpu::KernelInfo const&,
                    BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>> >
void HostDeviceFor (Gpu::KernelInfo const&,
                    BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void HostDeviceFor (Gpu::KernelInfo const&,
                    BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                    BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral_v<T1>>,
          typename M2=std::enable_if_t<std::is_integral_v<T2>>,
          typename M3=std::enable_if_t<std::is_integral_v<T3>> >
void HostDeviceFor (Gpu::KernelInfo const&,
                    BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                    BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    amrex::ignore_unused(MT);
    For(box1,ncomp1,std::forward<L1>(f1),
        box2,ncomp2,std::forward<L2>(f2),
        box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral_v<T>> >
AMREX_ATTRIBUTE_FLATTEN_FOR
void ParallelForRNG (T n, L const& f) noexcept
{
    for (T i = 0; i < n; ++i) {
        f(i,RandomEngine{});
    }
}

namespace detail {

template <int idim, typename L, int dim>
AMREX_FORCE_INLINE
void ParallelForRNG_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv) noexcept
{
    if constexpr (idim == 1) {
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_engine(f,iv,RandomEngine{});
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_engine(f,iv,RandomEngine{});
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_engine(f,iv,RandomEngine{});
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            ParallelForRNG_impND<idim-1>(f, lo, hi, iv);
        }
    }
}

template <int idim, typename L, typename T, int dim>
AMREX_FORCE_INLINE
void ParallelForRNG_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv, T n) noexcept
{
    if constexpr (idim == 1) {
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_engine(f,iv,n,RandomEngine{});
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_engine(f,iv,n,RandomEngine{});
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_engine(f,iv,n,RandomEngine{});
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            ParallelForRNG_impND<idim-1>(f, lo, hi, iv, n);
        }
    }
}

}

template <typename L, int dim>
AMREX_ATTRIBUTE_FLATTEN_FOR
void ParallelForRNG (BoxND<dim> const& box, L const& f) noexcept
{
    const auto lo = amrex::lbound_iv(box);
    const auto hi = amrex::ubound_iv(box);
    IntVectND<dim> iv;
    detail::ParallelForRNG_impND<dim>(f, lo, hi, iv);
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral_v<T>> >
AMREX_ATTRIBUTE_FLATTEN_FOR
void ParallelForRNG (BoxND<dim> const& box, T ncomp, L const& f) noexcept
{
    const auto lo = amrex::lbound_iv(box);
    const auto hi = amrex::ubound_iv(box);
    IntVectND<dim> iv;
    for (T n = 0; n < ncomp; ++n) {
        detail::ParallelForRNG_impND<dim>(f, lo, hi, iv, n);
    }
}

template <typename L>
void single_task (L&& f) noexcept
{
    std::forward<L>(f)();
}

}

#endif
