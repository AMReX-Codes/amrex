#ifndef AMREX_STENCIL_FOR_H_
#define AMREX_STENCIL_FOR_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_GpuTensorCopy.H>
#include <AMReX_IntVect.H>
#include <AMReX_TypeTraits.H>

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

namespace amrex::Gpu {

enum class TensorCopyPolicy {
    Auto,
    Always,
    Never
};

struct StencilInfo {
    IntVect tile = IntVect(AMREX_D_DECL(32, 8, 4));
    IntVect halo = IntVect(1);
    IntVect shared_padding = IntVect(0);
    TensorCopyPolicy tensor_copy = TensorCopyPolicy::Auto;
    bool require_full_tile = false;
    bool stage_one_component = true;

    StencilInfo& setTile (IntVect const& v) noexcept {
        tile = v;
        return *this;
    }

    StencilInfo& setHalo (IntVect const& v) noexcept {
        halo = v;
        return *this;
    }

    StencilInfo& setSharedPadding (IntVect const& v) noexcept {
        shared_padding = v;
        return *this;
    }

    StencilInfo& setTensorCopyPolicy (TensorCopyPolicy v) noexcept {
        tensor_copy = v;
        return *this;
    }

    StencilInfo& setRequireFullTile (bool v) noexcept {
        require_full_tile = v;
        return *this;
    }
};

template <typename T>
class TileView
{
public:
    AMREX_GPU_HOST_DEVICE
    constexpr TileView () noexcept = default;

    template <class U=T, std::enable_if_t<std::is_const_v<U>, int> = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr TileView (TileView<std::remove_const_t<T>> const& rhs) noexcept
        : m_p(rhs.m_p),
          m_stride_y(rhs.m_stride_y),
          m_stride_z(rhs.m_stride_z),
          m_stride_n(rhs.m_stride_n),
          m_ncomp(rhs.m_ncomp),
          m_tile_box(rhs.m_tile_box),
          m_valid_box(rhs.m_valid_box)
    {}

    AMREX_GPU_HOST_DEVICE
    constexpr TileView (T* p,
                        Long stride_y,
                        Long stride_z,
                        Long stride_n,
                        int ncomp,
                        Box const& tile_box,
                        Box const& valid_box) noexcept
        : m_p(p),
          m_stride_y(stride_y),
          m_stride_z(stride_z),
          m_stride_n(stride_n),
          m_ncomp(ncomp),
          m_tile_box(tile_box),
          m_valid_box(valid_box)
    {}

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T& operator() (int i, int j, int k) const noexcept {
        return this->operator()(i, j, k, 0);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T& operator() (int i, int j, int k, int n) const noexcept {
        AMREX_ASSERT(contains(i, j, k));
        AMREX_ASSERT(n >= 0 && n < m_ncomp);
        auto const lo = m_valid_box.smallEnd().dim3();
        return m_p[(i - lo.x)
                   + (j - lo.y) * m_stride_y
                   + (k - lo.z) * m_stride_z
                   + static_cast<Long>(n) * m_stride_n];
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool contains (int i, int j, int k) const noexcept {
        return m_valid_box.contains(i, j, k);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box tileBox () const noexcept {
        return m_tile_box;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box validBox () const noexcept {
        return m_valid_box;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int nComp () const noexcept {
        return m_ncomp;
    }

private:
    template <typename> friend class TileView;

    T* m_p = nullptr;
    Long m_stride_y = 1;
    Long m_stride_z = 1;
    Long m_stride_n = 1;
    int m_ncomp = 1;
    Box m_tile_box;
    Box m_valid_box;
};

} // namespace amrex::Gpu

namespace amrex::detail {

template <typename T>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Box array4_box_like (Array4<T> const& a, Box const& like) noexcept
{
    IntVect const lo(AMREX_D_DECL(a.begin.vect[0], a.begin.vect[1], a.begin.vect[2]));
    IntVect const hi(AMREX_D_DECL(a.end.vect[0]-1, a.end.vect[1]-1, a.end.vect[2]-1));
    return Box(lo, hi, like.ixType());
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
IntVect sanitize_stencil_tile (Box const& bx, Gpu::StencilInfo const& info) noexcept
{
    IntVect r = info.tile;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        r[d] = amrex::max(1, amrex::min(r[d], bx.length(d)));
    }
    return r;
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<int,3> tile_counts (Box const& bx, IntVect const& tile) noexcept
{
    auto const len = bx.length3d();
    Dim3 const t = tile.dim3(1);
    return {{(len[0] + t.x - 1) / t.x,
             (len[1] + t.y - 1) / t.y,
             (len[2] + t.z - 1) / t.z}};
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool tiles_are_uniform (Box const& bx, IntVect const& tile) noexcept
{
    auto const len = bx.length3d();
    Dim3 const t = tile.dim3(1);
    return (len[0] % t.x == 0) && (len[1] % t.y == 0) && (len[2] % t.z == 0);
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Box make_stencil_tile_box (Box const& bx,
                           IntVect const& tile,
                           int tx,
                           int ty,
                           int tz) noexcept
{
    IntVect small = bx.smallEnd() + IntVect(AMREX_D_DECL(tx*tile[0], ty*tile[1], tz*tile[2]));
    IntVect big = amrex::min(small + tile - IntVect(1), bx.bigEnd());
    return Box(small, big, bx.ixType());
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Box make_stencil_tile_box (Box const& bx,
                           IntVect const& tile,
                           int tile_id,
                           GpuArray<int,3> const& counts) noexcept
{
    int const tx = tile_id % counts[0];
    int const ty = (tile_id / counts[0]) % counts[1];
    int const tz = tile_id / (counts[0] * counts[1]);
    return make_stencil_tile_box(bx, tile, tx, ty, tz);
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Long num_tile_elements (IntVect const& extents) noexcept
{
    Dim3 const e = extents.dim3(1);
    return static_cast<Long>(e.x) * static_cast<Long>(e.y) * static_cast<Long>(e.z);
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Long num_tile_elements (Box const& bx) noexcept
{
    return num_tile_elements(bx.length());
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
GpuArray<Long,3> tile_view_strides (Box const& valid_box,
                                    IntVect const& padding,
                                    int ncomp = 1) noexcept
{
    amrex::ignore_unused(ncomp);
    Dim3 const ext = valid_box.length().dim3(1);
    Dim3 const pad = padding.dim3(0);
    Long const stride_y = static_cast<Long>(ext.x + pad.x);
    Long const stride_z = stride_y * static_cast<Long>(ext.y + pad.y);
    Long const stride_n = stride_z * static_cast<Long>(ext.z + pad.z);
    return {{stride_y, stride_z, stride_n}};
}

template <typename T>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Gpu::TileView<T> make_tile_view (T* p,
                                 Box const& tile_box,
                                 Box const& valid_box,
                                 IntVect const& padding,
                                 int ncomp = 1) noexcept
{
    auto const strides = tile_view_strides(valid_box, padding, ncomp);
    return Gpu::TileView<T>(p, strides[0], strides[1], strides[2], ncomp, tile_box, valid_box);
}

template <typename T>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Gpu::TileView<T> make_tile_view (T* p,
                                 Box const& tile_box,
                                 Box const& valid_box,
                                 Long stride_y,
                                 Long stride_z,
                                 Long stride_n,
                                 int ncomp = 1) noexcept
{
    return Gpu::TileView<T>(p, stride_y, stride_z, stride_n, ncomp, tile_box, valid_box);
}

template <typename T>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Gpu::TileView<T const> make_global_alias_tile_view (Array4<T const> const& src,
                                                    Box const& tile_box,
                                                    Box const& valid_box) noexcept
{
    auto const lo = valid_box.smallEnd().dim3();
    return make_tile_view(src.ptr(lo.x, lo.y, lo.z), tile_box, valid_box,
                          src.stride.a[0], src.stride.a[1], src.stride.a[2], 1);
}

template <typename T>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::size_t shared_tile_bytes (IntVect const& valid_extents,
                               IntVect const& padding,
                               int ncomp = 1) noexcept
{
    Dim3 const ext = valid_extents.dim3(1);
    Dim3 const pad = padding.dim3(0);
    auto const nelems = static_cast<std::size_t>(ext.x + pad.x)
        * static_cast<std::size_t>(ext.y + pad.y)
        * static_cast<std::size_t>(ext.z + pad.z)
        * static_cast<std::size_t>(ncomp);
    return nelems * sizeof(T);
}

template <typename T>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::size_t shared_tile_bytes (Box const& valid_box,
                               IntVect const& padding,
                               int ncomp = 1) noexcept
{
    return shared_tile_bytes<T>(valid_box.length(), padding, ncomp);
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool boxes_equal_bounds (Box const& lhs, Box const& rhs) noexcept
{
    return lhs.smallEnd() == rhs.smallEnd() && lhs.bigEnd() == rhs.bigEnd();
}

template <typename F, typename Tile>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void call_stencil_body (F const& f, Tile const& tile, int i, int j, int k) noexcept
{
    f(tile, i, j, k);
}

template <int MT, typename T, typename F>
void StencilFor_doit (Box const& bx,
                      Gpu::StencilInfo const& info,
                      Array4<T const> const& src,
                      F&& f) noexcept;

} // namespace amrex::detail

#ifdef AMREX_USE_GPU
#include <AMReX_StencilForG.H>
#else
#include <AMReX_StencilForC.H>
#endif

namespace amrex {

template <int MT = AMREX_GPU_MAX_THREADS, typename T, typename F>
void StencilFor (Box const& bx,
                 Gpu::StencilInfo const& info,
                 Array4<T const> const& src,
                 F&& f) noexcept
{
    detail::StencilFor_doit<MT>(bx, info, src, std::forward<F>(f));
}

template <int MT = AMREX_GPU_MAX_THREADS, typename T, typename F>
void StencilFor (Box const& bx,
                 int ncomp,
                 Gpu::StencilInfo const& info,
                 Array4<T const> const& src,
                 F&& f) noexcept
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncomp >= 0, "StencilFor: ncomp must be non-negative");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(src.nComp() >= ncomp,
                                     "StencilFor: staged source does not have enough components");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(info.stage_one_component,
                                     "StencilFor: multi-component staging is not implemented in this prototype");

    for (int n = 0; n < ncomp; ++n) {
        auto const src_n = Array4<T const>(src, n, 1);
        StencilFor<MT>(bx, info, src_n,
        [=] AMREX_GPU_DEVICE (auto const& tile, int i, int j, int k) noexcept
        {
            f(tile, i, j, k, n);
        });
    }
}

} // namespace amrex

#endif
