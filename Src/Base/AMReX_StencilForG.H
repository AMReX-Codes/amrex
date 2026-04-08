#ifndef AMREX_STENCIL_FOR_G_H_
#define AMREX_STENCIL_FOR_G_H_
#include <AMReX_Config.H>

#ifdef AMREX_USE_GPU

#include <AMReX_GpuError.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuLaunchFunctsG.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Loop.H>

#if AMREX_USE_CUDA_BLOCK_LOAD_TO_SHARED
#include <cub/block/block_load_to_shared.cuh>
#include <cuda/std/span>
#endif

#include <algorithm>
#include <type_traits>

namespace amrex::detail {

template <typename T, typename F>
void StencilFor_direct_alias_fallback (Box const& bx,
                                       IntVect const& tile,
                                       Gpu::StencilInfo const& info,
                                       Array4<T const> const& src,
                                       F const& f) noexcept
{
    Box const src_box = array4_box_like(src, bx);
    auto const counts = tile_counts(bx, tile);

    for (int tz = 0; tz < counts[2]; ++tz) {
        for (int ty = 0; ty < counts[1]; ++ty) {
            for (int tx = 0; tx < counts[0]; ++tx) {
                Box const tile_box = make_stencil_tile_box(bx, tile, tx, ty, tz);
                Box const staged_box = amrex::grow(tile_box, info.halo);
                Box const valid_box = staged_box & src_box;

                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!info.require_full_tile || boxes_equal_bounds(valid_box, staged_box),
                                                 "StencilFor: full staged tile is not available");

                auto const tile_view = make_global_alias_tile_view(src, tile_box, valid_box);
                amrex::ParallelFor(tile_box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    call_stencil_body(f, tile_view, i, j, k);
                });
            }
        }
    }
}

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

template <typename T>
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
char* align_shared_ptr (char* p) noexcept
{
    auto const a = static_cast<std::uintptr_t>(Gpu::detail::tensorCopyAlignment());
    auto const v = reinterpret_cast<std::uintptr_t>(p);
    return reinterpret_cast<char*>((v + (a - 1u)) & ~(a - 1u));
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void stencil_copy_manual_to_shared (Array4<T const> const& src,
                                    Box const& valid_box,
                                    T* shared_ptr,
                                    Long stride_y,
                                    Long stride_z) noexcept
{
    auto const ext = valid_box.length3d();
    Long const npts = static_cast<Long>(ext[0]) * static_cast<Long>(ext[1]) * static_cast<Long>(ext[2]);
    Dim3 const vlo = valid_box.smallEnd().dim3();

    for (Long idx = static_cast<Long>(threadIdx.x); idx < npts; idx += static_cast<Long>(blockDim.x)) {
        Long tmp = idx;
        int const ii = static_cast<int>(tmp % ext[0]);
        tmp /= ext[0];
        int const jj = static_cast<int>(tmp % ext[1]);
        int const kk = static_cast<int>(tmp / ext[1]);

        int const i = vlo.x + ii;
        int const j = vlo.y + jj;
        int const k = vlo.z + kk;
        shared_ptr[static_cast<Long>(ii) + static_cast<Long>(jj) * stride_y + static_cast<Long>(kk) * stride_z]
            = src(i, j, k, 0);
    }
}

#if AMREX_USE_CUDA_BLOCK_LOAD_TO_SHARED
template <int MT, typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void stencil_copy_tensor_to_shared (Array4<T const> const& src,
                                    Box const& valid_box,
                                    T* shared_ptr,
                                    Long stride_y,
                                    Long stride_z) noexcept
{
    using Loader = cub::detail::BlockLoadToShared<MT>;
    __shared__ typename Loader::TempStorage load_temp_storage;

    Loader loader(load_temp_storage);
    auto const ext = valid_box.length3d();
    Dim3 const vlo = valid_box.smallEnd().dim3();

    for (int kk = 0; kk < ext[2]; ++kk) {
        for (int jj = 0; jj < ext[1]; ++jj) {
            T* const dst = shared_ptr + static_cast<Long>(jj) * stride_y + static_cast<Long>(kk) * stride_z;
            T const* const src_row = src.ptr(vlo.x, vlo.y + jj, vlo.z + kk);
            auto const dst_span = cuda::std::span<T>(dst, ext[0]);
            auto const src_span = cuda::std::span<T const>(src_row, ext[0]);
            loader.template CopyAsync<T, Gpu::detail::tensorCopyAlignment()>(dst_span, src_span);
        }
    }

    auto token = loader.Commit();
    loader.Wait(static_cast<typename Loader::CommitToken&&>(token));
}
#endif

template <int MT, bool UseTensorCopy, typename T, typename F>
void StencilFor_shared_launch (Box const& bx,
                               IntVect const& tile,
                               Gpu::StencilInfo const& info,
                               Array4<T const> const& src,
                               F const& f,
                               std::size_t shared_mem_bytes) noexcept
{
    auto const counts = tile_counts(bx, tile);
    int const ntile_blocks = counts[0] * counts[1] * counts[2];

    amrex::launch_global<MT><<<ntile_blocks, MT, shared_mem_bytes, Gpu::gpuStream()>>>(
    [=] AMREX_GPU_DEVICE () noexcept
    {
        int const tile_id = static_cast<int>(blockIdx.x);
        if (tile_id >= ntile_blocks) {
            return;
        }

        Box const src_box = array4_box_like(src, bx);
        Box const tile_box = make_stencil_tile_box(bx, tile, tile_id, counts);
        Box const valid_box = amrex::grow(tile_box, info.halo) & src_box;
        auto const strides = tile_view_strides(valid_box, info.shared_padding);

        Gpu::SharedMemory<char> smem;
        char* const raw_shared = smem.dataPtr();
        T* const shared_ptr = reinterpret_cast<T*>(align_shared_ptr<T>(raw_shared));

        if constexpr (UseTensorCopy) {
#if AMREX_USE_CUDA_BLOCK_LOAD_TO_SHARED
            stencil_copy_tensor_to_shared<MT>(src, valid_box, shared_ptr, strides[0], strides[1]);
#else
            stencil_copy_manual_to_shared(src, valid_box, shared_ptr, strides[0], strides[1]);
            __syncthreads();
#endif
        } else {
            stencil_copy_manual_to_shared(src, valid_box, shared_ptr, strides[0], strides[1]);
            __syncthreads();
        }

        auto const tile_view = make_tile_view(shared_ptr, tile_box, valid_box, strides[0], strides[1], strides[2], 1);
        auto const ext = tile_box.length3d();
        Long const npts = static_cast<Long>(ext[0]) * static_cast<Long>(ext[1]) * static_cast<Long>(ext[2]);
        Dim3 const tlo = tile_box.smallEnd().dim3();

        for (Long idx = static_cast<Long>(threadIdx.x); idx < npts; idx += static_cast<Long>(blockDim.x)) {
            Long tmp = idx;
            int const ii = static_cast<int>(tmp % ext[0]);
            tmp /= ext[0];
            int const jj = static_cast<int>(tmp % ext[1]);
            int const kk = static_cast<int>(tmp / ext[1]);
            call_stencil_body(f, tile_view, tlo.x + ii, tlo.y + jj, tlo.z + kk);
        }
    });

    AMREX_GPU_ERROR_CHECK();
}

template <int MT, typename T>
[[nodiscard]] inline bool can_use_native_tensor_copy (Box const& bx,
                                                      IntVect const& tile,
                                                      Gpu::StencilInfo const& info,
                                                      Array4<T const> const& src,
                                                      std::size_t shared_mem_bytes) noexcept
{
    if (info.tensor_copy == Gpu::TensorCopyPolicy::Never) {
        return false;
    }

#if AMREX_USE_CUDA_BLOCK_LOAD_TO_SHARED
    if constexpr (MT % Gpu::Device::warp_size != 0) {
        return false;
    }

    if (!Gpu::detail::queryTensorCopyCaps().available) {
        return false;
    }

    if (!tiles_are_uniform(bx, tile)) {
        return false;
    }

    Box const src_box = array4_box_like(src, bx);
    if (!src_box.contains(amrex::grow(bx, info.halo))) {
        return false;
    }

    if (shared_mem_bytes > Gpu::Device::sharedMemPerBlock()) {
        return false;
    }

    Box const full_stage = amrex::grow(make_stencil_tile_box(bx, tile, 0, 0, 0), info.halo);
    auto const stage_lo = full_stage.smallEnd().dim3();
    T const* const stage_ptr = src.ptr(stage_lo.x, stage_lo.y, stage_lo.z);
    Long const row_bytes = static_cast<Long>(full_stage.length(0)) * sizeof(T);

    return Gpu::detail::canUseTensorCopy(stage_ptr, src_box, tile, info.halo, sizeof(T))
        && (row_bytes % static_cast<Long>(Gpu::detail::tensorCopyAlignment()) == 0)
        && Gpu::detail::isTensorCopyAligned(stage_ptr)
        && ((src.stride.a[0] * sizeof(T)) % static_cast<Long>(Gpu::detail::tensorCopyAlignment()) == 0)
        && ((src.stride.a[1] * sizeof(T)) % static_cast<Long>(Gpu::detail::tensorCopyAlignment()) == 0);
#else
    amrex::ignore_unused(bx, tile, info, src, shared_mem_bytes);
    return false;
#endif
}

#endif

template <int MT, typename T, typename F>
void StencilFor_doit (Box const& bx,
                      Gpu::StencilInfo const& info,
                      Array4<T const> const& src,
                      F&& f) noexcept
{
    if (bx.isEmpty()) {
        return;
    }

    IntVect const tile = sanitize_stencil_tile(bx, info);
    Box const src_box = array4_box_like(src, bx);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(src_box.contains(bx),
                                     "StencilFor: source Array4 must contain the compute box");

#if defined(AMREX_USE_SYCL)
    StencilFor_direct_alias_fallback(bx, tile, info, src, f);
#elif defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(MT <= Gpu::Device::maxThreadsPerBlock(),
                                     "StencilFor: MT exceeds the device thread-block limit");

    IntVect const shared_extents = tile + info.halo + info.halo;
    std::size_t const shared_mem_bytes = shared_tile_bytes<T>(shared_extents, info.shared_padding)
        + Gpu::detail::tensorCopyAlignment();

    if (shared_mem_bytes > Gpu::Device::sharedMemPerBlock()) {
        StencilFor_direct_alias_fallback(bx, tile, info, src, f);
        return;
    }

    if (can_use_native_tensor_copy<MT, T>(bx, tile, info, src, shared_mem_bytes)) {
        StencilFor_shared_launch<MT, true>(bx, tile, info, src, f, shared_mem_bytes);
    } else {
        StencilFor_shared_launch<MT, false>(bx, tile, info, src, f, shared_mem_bytes);
    }
#else
    StencilFor_direct_alias_fallback(bx, tile, info, src, f);
#endif
}

} // namespace amrex::detail

#endif
#endif
