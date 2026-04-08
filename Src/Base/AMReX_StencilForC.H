#ifndef AMREX_STENCIL_FOR_C_H_
#define AMREX_STENCIL_FOR_C_H_
#include <AMReX_Config.H>

#ifndef AMREX_USE_GPU

#include <AMReX_Loop.H>
#include <AMReX_Vector.H>

namespace amrex::detail {

template <int MT, typename T, typename F>
void StencilFor_doit (Box const& bx,
                      Gpu::StencilInfo const& info,
                      Array4<T const> const& src,
                      F&& f) noexcept
{
    amrex::ignore_unused(MT);

    if (bx.isEmpty()) {
        return;
    }

    IntVect const tile = sanitize_stencil_tile(bx, info);
    Box const src_box = array4_box_like(src, bx);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(src_box.contains(bx),
                                     "StencilFor: source Array4 must contain the compute box");

    auto const counts = tile_counts(bx, tile);

    for (int tz = 0; tz < counts[2]; ++tz) {
        for (int ty = 0; ty < counts[1]; ++ty) {
            for (int tx = 0; tx < counts[0]; ++tx) {
                Box const tile_box = make_stencil_tile_box(bx, tile, tx, ty, tz);
                Box const staged_box = amrex::grow(tile_box, info.halo);
                Box const valid_box = staged_box & src_box;

                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!info.require_full_tile || boxes_equal_bounds(valid_box, staged_box),
                                                 "StencilFor: full staged tile is not available");

                Vector<T> tile_storage(static_cast<std::size_t>(shared_tile_bytes<T>(valid_box, info.shared_padding))
                                       / sizeof(T));
                auto stage_view = make_tile_view(tile_storage.data(), tile_box, valid_box, info.shared_padding);

                amrex::LoopOnCpu(valid_box, [&] (int i, int j, int k) noexcept
                {
                    stage_view(i, j, k) = src(i, j, k, 0);
                });

                amrex::LoopOnCpu(tile_box, [&] (int i, int j, int k) noexcept
                {
                    call_stencil_body(f, stage_view, i, j, k);
                });
            }
        }
    }
}

} // namespace amrex::detail

#endif
#endif
