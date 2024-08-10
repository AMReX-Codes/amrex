#ifndef AMREX_LOOP_H_
#define AMREX_LOOP_H_
#include <AMReX_Config.H>

#include <AMReX_Box.H>
#include <AMReX_Extension.H>

namespace amrex {

namespace loop_detail {

    // call_f_intvect_inner

    template <typename F, std::size_t...Ns, class...Args>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<1> iv, Args...args)
        noexcept -> decltype(f(0, 0, 0, args...))
    {
        f(iv[0], 0, 0, args...);
    }

    template <typename F, std::size_t...Ns, class...Args>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<2> iv, Args...args)
        noexcept -> decltype(f(0, 0, 0, args...))
    {
        f(iv[0], iv[1], 0, args...);
    }

    template <typename F, int dim, std::size_t...Ns, class...Args>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<dim> iv, Args...args)
        noexcept -> decltype(f(iv, args...))
    {
        f(iv, args...);
    }

    template <typename F, int dim, std::size_t...Ns, class...Args>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<dim> iv, Args...args)
        noexcept -> decltype(f(iv[Ns]..., args...))
    {
        f(iv[Ns]..., args...);
    }

    // call_f_intvect

    template <typename F, int dim>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect (F const& f, IntVectND<dim> iv)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv);
    }

    // call_f_intvect_ncomp

    template <typename F, int dim>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_ncomp (F const& f, IntVectND<dim> iv, int n)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, 0))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n);
    }

    // call_f_intvect_inner_cpu

    template <typename F, std::size_t...Ns, class...Args>
    AMREX_FORCE_INLINE
    auto call_f_intvect_inner_cpu (std::index_sequence<Ns...>, F const& f, IntVectND<1> iv, Args...args)
        noexcept -> decltype(f(0, 0, 0, args...))
    {
        f(iv[0], 0, 0, args...);
    }

    template <typename F, std::size_t...Ns, class...Args>
    AMREX_FORCE_INLINE
    auto call_f_intvect_inner_cpu (std::index_sequence<Ns...>, F const& f, IntVectND<2> iv, Args...args)
        noexcept -> decltype(f(0, 0, 0, args...))
    {
        f(iv[0], iv[1], 0, args...);
    }

    template <typename F, int dim, std::size_t...Ns, class...Args>
    AMREX_FORCE_INLINE
    auto call_f_intvect_inner_cpu (std::index_sequence<Ns...>, F const& f, IntVectND<dim> iv, Args...args)
        noexcept -> decltype(f(iv, args...))
    {
        f(iv, args...);
    }

    template <typename F, int dim, std::size_t...Ns, class...Args>
    AMREX_FORCE_INLINE
    auto call_f_intvect_inner_cpu (std::index_sequence<Ns...>, F const& f, IntVectND<dim> iv, Args...args)
        noexcept -> decltype(f(iv[Ns]..., args...))
    {
        f(iv[Ns]..., args...);
    }

    // call_f_intvect_cpu

    template <typename F, int dim>
    AMREX_FORCE_INLINE
    auto call_f_intvect_cpu (F const& f, IntVectND<dim> iv)
        noexcept -> decltype(call_f_intvect_inner_cpu(std::make_index_sequence<dim>(), f, iv))
    {
        call_f_intvect_inner_cpu(std::make_index_sequence<dim>(), f, iv);
    }

    // call_f_intvect_ncomp_cpu

    template <typename F, int dim>
    AMREX_FORCE_INLINE
    auto call_f_intvect_ncomp_cpu (F const& f, IntVectND<dim> iv, int n)
        noexcept -> decltype(call_f_intvect_inner_cpu(std::make_index_sequence<dim>(), f, iv, 0))
    {
        call_f_intvect_inner_cpu(std::make_index_sequence<dim>(), f, iv, n);
    }
}

template <class F>
AMREX_GPU_HOST_DEVICE
AMREX_ATTRIBUTE_FLATTEN_FOR
void Loop (Dim3 lo, Dim3 hi, F const& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_GPU_HOST_DEVICE
AMREX_ATTRIBUTE_FLATTEN_FOR
void Loop (Dim3 lo, Dim3 hi, int ncomp, F const& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_GPU_HOST_DEVICE
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopConcurrent (Dim3 lo, Dim3 hi, F const& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_GPU_HOST_DEVICE
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopConcurrent (Dim3 lo, Dim3 hi, int ncomp, F const& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

namespace loop_detail {

template <int idim, typename L, int dim>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Loop_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv) noexcept
{
    if constexpr (idim == 1) {
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect(f,iv);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect(f,iv);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect(f,iv);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            Loop_impND<idim-1>(f, lo, hi, iv);
        }
    }
}

}

template <class F, int dim>
AMREX_GPU_HOST_DEVICE
AMREX_ATTRIBUTE_FLATTEN_FOR
void Loop (BoxND<dim> const& bx, F const& f) noexcept
{
    const auto lo = amrex::lbound_iv(bx);
    const auto hi = amrex::ubound_iv(bx);
    IntVectND<dim> iv;
    loop_detail::Loop_impND<dim>(f, lo, hi, iv);
}

namespace loop_detail {

template <int idim, typename L, int dim>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Loop_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv, int n) noexcept
{
    if constexpr (idim == 1) {
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp(f,iv,n);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp(f,iv,n);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp(f,iv,n);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            Loop_impND<idim-1>(f, lo, hi, iv, n);
        }
    }
}

}

template <class F, int dim>
AMREX_GPU_HOST_DEVICE
AMREX_ATTRIBUTE_FLATTEN_FOR
void Loop (BoxND<dim> const& bx, int ncomp, F const& f) noexcept
{
    const auto lo = amrex::lbound_iv(bx);
    const auto hi = amrex::ubound_iv(bx);
    IntVectND<dim> iv;
    for (int n = 0; n < ncomp; ++n) {
        loop_detail::Loop_impND<dim>(f, lo, hi, iv, n);
    }
}

namespace loop_detail {

template <int idim, typename L, int dim>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void LoopConcurrent_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv) noexcept
{
    if constexpr (idim == 1) {
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect(f,iv);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect(f,iv);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect(f,iv);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            LoopConcurrent_impND<idim-1>(f, lo, hi, iv);
        }
    }
}

}

template <class F, int dim>
AMREX_GPU_HOST_DEVICE
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopConcurrent (BoxND<dim> const& bx, F const& f) noexcept
{
    const auto lo = amrex::lbound_iv(bx);
    const auto hi = amrex::ubound_iv(bx);
    IntVectND<dim> iv;
    loop_detail::LoopConcurrent_impND<dim>(f, lo, hi, iv);
}

namespace loop_detail {

template <int idim, typename L, int dim>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void LoopConcurrent_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv, int n) noexcept
{
    if constexpr (idim == 1) {
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp(f,iv,n);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp(f,iv,n);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp(f,iv,n);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            LoopConcurrent_impND<idim-1>(f, lo, hi, iv, n);
        }
    }
}

}

template <class F, int dim>
AMREX_GPU_HOST_DEVICE
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopConcurrent (BoxND<dim> const& bx, int ncomp, F const& f) noexcept
{
    const auto lo = amrex::lbound_iv(bx);
    const auto hi = amrex::ubound_iv(bx);
    IntVectND<dim> iv;
    for (int n = 0; n < ncomp; ++n) {
        loop_detail::LoopConcurrent_impND<dim>(f, lo, hi, iv, n);
    }
}

// The functions above are __host__ __device__ functions.  If f is not a
// __host__ __device__ function, we will get warning about calling __host__
// function from a __host__ __device__ function.  This is ugly.  To get rid
// of the warning, we have to use the functions below for those situations.

template <class F>
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopOnCpu (Dim3 lo, Dim3 hi, F const& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopOnCpu (Dim3 lo, Dim3 hi, int ncomp, F const& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

template <class F>
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopConcurrentOnCpu (Dim3 lo, Dim3 hi, F const& f) noexcept
{
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k);
    }}}
}

template <class F>
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopConcurrentOnCpu (Dim3 lo, Dim3 hi, int ncomp, F const& f) noexcept
{
    for (int n = 0; n < ncomp; ++n) {
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        f(i,j,k,n);
    }}}}
}

namespace loop_detail {

template <int idim, typename L, int dim>
AMREX_FORCE_INLINE
void LoopOnCpu_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv) noexcept
{
    if constexpr (idim == 1) {
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_cpu(f,iv);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_cpu(f,iv);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_cpu(f,iv);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            LoopOnCpu_impND<idim-1>(f, lo, hi, iv);
        }
    }
}

}

template <class F, int dim>
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopOnCpu (BoxND<dim> const& bx, F const& f) noexcept
{
    const auto lo = amrex::lbound_iv(bx);
    const auto hi = amrex::ubound_iv(bx);
    IntVectND<dim> iv;
    loop_detail::LoopOnCpu_impND<dim>(f, lo, hi, iv);
}

namespace loop_detail {

template <int idim, typename L, int dim>
AMREX_FORCE_INLINE
void LoopOnCpu_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv, int n) noexcept
{
    if constexpr (idim == 1) {
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_cpu(f,iv,n);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_cpu(f,iv,n);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_cpu(f,iv,n);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            LoopOnCpu_impND<idim-1>(f, lo, hi, iv, n);
        }
    }
}

}

template <class F, int dim>
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopOnCpu (BoxND<dim> const& bx, int ncomp, F const& f) noexcept
{
    const auto lo = amrex::lbound_iv(bx);
    const auto hi = amrex::ubound_iv(bx);
    IntVectND<dim> iv;
    for (int n = 0; n < ncomp; ++n) {
        loop_detail::LoopOnCpu_impND<dim>(f, lo, hi, iv, n);
    }
}

namespace loop_detail {

template <int idim, typename L, int dim>
AMREX_FORCE_INLINE
void LoopConcurrentOnCpu_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv) noexcept
{
    if constexpr (idim == 1) {
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_cpu(f,iv);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_cpu(f,iv);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_cpu(f,iv);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            LoopConcurrentOnCpu_impND<idim-1>(f, lo, hi, iv);
        }
    }
}

}

template <class F, int dim>
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopConcurrentOnCpu (BoxND<dim> const& bx, F const& f) noexcept
{
    const auto lo = amrex::lbound_iv(bx);
    const auto hi = amrex::ubound_iv(bx);
    IntVectND<dim> iv;
    loop_detail::LoopConcurrentOnCpu_impND<dim>(f, lo, hi, iv);
}

namespace loop_detail {

template <int idim, typename L, int dim>
AMREX_FORCE_INLINE
void LoopConcurrentOnCpu_impND (L const& f, IntVectND<dim> const lo, IntVectND<dim> const hi, IntVectND<dim> iv, int n) noexcept
{
    if constexpr (idim == 1) {
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_cpu(f,iv,n);
        }
    } else if constexpr (idim == 2) {
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_cpu(f,iv,n);
        }}
    } else if constexpr (idim == 3) {
        for (int i2 = lo[2], h2 = hi[2]; i2 <= h2; ++i2) { iv[2] = i2;
        for (int i1 = lo[1], h1 = hi[1]; i1 <= h1; ++i1) { iv[1] = i1;
        AMREX_PRAGMA_SIMD
        for (int i0 = lo[0], h0 = hi[0]; i0 <= h0; ++i0) { iv[0] = i0;
            call_f_intvect_ncomp_cpu(f,iv,n);
        }}}
    } else {
        for (int id = lo[idim-1], hd = hi[idim-1]; id <= hd; ++id) { iv[idim-1] = id;
            LoopConcurrentOnCpu_impND<idim-1>(f, lo, hi, iv, n);
        }
    }
}

}

template <class F, int dim>
AMREX_ATTRIBUTE_FLATTEN_FOR
void LoopConcurrentOnCpu (BoxND<dim> const& bx, int ncomp, F const& f) noexcept
{
    const auto lo = amrex::lbound_iv(bx);
    const auto hi = amrex::ubound_iv(bx);
    IntVectND<dim> iv;
    for (int n = 0; n < ncomp; ++n) {
        loop_detail::LoopConcurrentOnCpu_impND<dim>(f, lo, hi, iv, n);
    }
}

// Implementation of "constexpr for" based on
// https://artificial-mind.net/blog/2020/10/31/constexpr-for
//
// Approximates what one would get from a compile-time
// unrolling of the loop
// for (int i = 0; i < N; ++i) {
//    f(i);
// }
//
// The mechanism is recursive: we evaluate f(i) at the current
// i and then call the for loop at i+1. f() is a lambda function
// that provides the body of the loop and takes only an integer
// i as its argument.

template<auto I, auto N, class F>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr void constexpr_for (F const& f)
{
    if constexpr (I < N) {
        f(std::integral_constant<decltype(I), I>());
        constexpr_for<I+1, N>(f);
    }
}

#include <AMReX_Loop.nolint.H>

}

#endif
