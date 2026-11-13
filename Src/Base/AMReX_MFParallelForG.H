#ifndef AMREX_MF_PARALLEL_FOR_G_H_
#define AMREX_MF_PARALLEL_FOR_G_H_
#include <AMReX_Config.H>

#ifdef AMREX_USE_GPU

#include <AMReX_Concepts.H>
#include <AMReX_GpuDevice.H>
#include <algorithm>
#include <cmath>
#include <limits>

/// \cond DOXYGEN_IGNORE
namespace amrex::detail {

// Copies the BoxIndexers of boxes, followed by the ints in extra, to the device.
inline
void build_par_for_boxes (char*& hp, BoxIndexer*& pboxes, Vector<Box> const& boxes,
                          Vector<int> const& extra = {})
{
    if (boxes.empty()) { return; }
    static_assert(sizeof(BoxIndexer) % alignof(int) == 0);
    const int nboxes = boxes.size();
    const std::size_t nbytes_boxes = nboxes*sizeof(BoxIndexer);
    const std::size_t nbytes = nbytes_boxes + extra.size()*sizeof(int);
    hp = (char*)The_Pinned_Arena()->alloc(nbytes);
    auto* hp_boxes = (BoxIndexer*)hp;
    for (int i = 0; i < nboxes; ++i) {
        new (hp_boxes+i) BoxIndexer(boxes[i]);
    }
    auto* hp_extra = (int*)(hp + nbytes_boxes);
    for (int i = 0; i < int(extra.size()); ++i) {
        hp_extra[i] = extra[i];
    }

    auto dp = (char*) The_Arena()->alloc(nbytes);
    Gpu::htod_memcpy_async(dp, hp, nbytes);
    pboxes = (BoxIndexer*)dp;
}

inline
void destroy_par_for_boxes (char* hp, char* dp)
{
    // The host-to-device copy in build_par_for_boxes() is asynchronous.  In a
    // NoSync region it may still be pending when the owning ParForInfo is
    // destroyed (e.g., when a temporary FabArray with a unique BoxArray dies).
    // Use stream-ordered frees so that the buffers are not reused too early.
    if (hp) {
        Gpu::freeAsync(The_Pinned_Arena(), hp);
        Gpu::freeAsync(The_Arena(), dp);
    }
}

namespace parfor_mf_detail {
    template <typename F>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int b, int i, int j, int k, int) noexcept
        -> decltype(f(0,0,0,0))
    {
        f(b,i,j,k);
    }

    template <typename F>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int b, int i, int j, int k, int ncomp) noexcept
        -> decltype(f(0,0,0,0,0))
    {
        for (int n = 0; n < ncomp; ++n) {
            f(b,i,j,k,n);
        }
    }
}

template <int MT, FabArrayType MF, typename F>
void
ParallelFor_doit (MF const& mf, IntVect const& nghost, int ncomp, IntVect const&, bool, F const& f)
{
    const auto& index_array = mf.IndexArray();
    const int nboxes = index_array.size();

    if (nboxes == 0) {
        return;
    } else if (nboxes == 1) {
        Box const& b = amrex::grow(mf.box(index_array[0]), nghost);
        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            parfor_mf_detail::call_f(f, 0, i, j, k, ncomp);
        });
    } else {
        auto const& parforinfo = mf.getParForInfo(nghost);
        auto nblocks_per_box = parforinfo.getNBlocksPerBox(MT);
        AMREX_ASSERT(Long(nblocks_per_box)*Long(nboxes) < Long(std::numeric_limits<int>::max()));
        const int nblocks = nblocks_per_box * nboxes;
        const BoxIndexer* dp_boxes = parforinfo.getBoxes();
        // Box and block within the box from the block index, without a division.
        const Math::FastDivmodU32 fdm_box(nblocks_per_box);

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

        amrex::launch_global<MT>
            <<<nblocks, MT, 0, Gpu::gpuStream()>>>
            ([=] AMREX_GPU_DEVICE () noexcept
             {
                 std::uint32_t ibox, iblock;
                 fdm_box(ibox, iblock, blockIdx.x);
                 auto icell = BoxIndexer::UInt(iblock)*MT + threadIdx.x;

#elif defined(AMREX_USE_SYCL)

        amrex::launch<MT>(nblocks, Gpu::gpuStream(),
             [=] AMREX_GPU_DEVICE (sycl::nd_item<1> const& item) noexcept
             {
                 auto blockIdxx = std::uint32_t(item.get_group_linear_id());
                 int threadIdxx = item.get_local_linear_id();
                 std::uint32_t ibox, iblock;
                 fdm_box(ibox, iblock, blockIdxx);
                 auto icell = BoxIndexer::UInt(iblock)*MT + threadIdxx;
#endif
                 BoxIndexer const& indexer = dp_boxes[ibox];
                 if (icell < indexer.numPts()) {
                     auto [i, j, k] = indexer(icell);
                     parfor_mf_detail::call_f(f, int(ibox), i, j, k, ncomp);
                 }
             });
    }
    AMREX_GPU_ERROR_CHECK();
}

template <FabArrayType MF, typename F>
void
ParallelFor_doit (MF const& mf, IntVect const& nghost, int ncomp, IntVect const& ts, bool dynamic, F&& f)
{
#ifdef AMREX_USE_CUDA
    constexpr int MT = 128;
#else
    constexpr int MT = AMREX_GPU_MAX_THREADS;
#endif
    ParallelFor_doit<MT>(mf, nghost, ncomp, ts, dynamic, std::forward<F>(f));
}

template <FabArrayType MF, typename F>
void
ParallelFor_doit (MF const& mf, IntVect const& nghost, IntVect const& ts, bool dynamic, F&& f)
{
#ifdef AMREX_USE_CUDA
    constexpr int MT = 128;
#else
    constexpr int MT = AMREX_GPU_MAX_THREADS;
#endif
    ParallelFor_doit<MT>(mf, nghost, 1, ts, dynamic, std::forward<F>(f));
}

// One launch over the valid points of all boxes whose index equals offset
// modulo stride.
template <int MT, FabArrayType MF, typename F>
void
ParallelForStrided_doit (MF const& mf, IntVect const& stride, IntVect const& offset, F const& f)
{
    const int nboxes = mf.IndexArray().size();
    if (nboxes == 0) { return; }

    auto const& parforinfo = mf.getParForInfo(stride, offset);
    const auto nblocks_per_box = parforinfo.getNBlocksPerBox(MT);
    if (nblocks_per_box == 0) { return; }
    AMREX_ASSERT(Long(nblocks_per_box)*Long(nboxes) < Long(std::numeric_limits<int>::max()));
    const int nblocks = nblocks_per_box * nboxes;
    const BoxIndexer* dp_boxes = parforinfo.getBoxes();
    const Math::FastDivmodU32 fdm_box(nblocks_per_box);

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    amrex::launch_global<MT>
        <<<nblocks, MT, 0, Gpu::gpuStream()>>>
        ([=] AMREX_GPU_DEVICE () noexcept
         {
             std::uint32_t ibox, iblock;
             fdm_box(ibox, iblock, blockIdx.x);
             auto icell = BoxIndexer::UInt(iblock)*MT + threadIdx.x;

#elif defined(AMREX_USE_SYCL)

    amrex::launch<MT>(nblocks, Gpu::gpuStream(),
         [=] AMREX_GPU_DEVICE (sycl::nd_item<1> const& item) noexcept
         {
             auto blockIdxx = std::uint32_t(item.get_group_linear_id());
             int threadIdxx = item.get_local_linear_id();
             std::uint32_t ibox, iblock;
             fdm_box(ibox, iblock, blockIdxx);
             auto icell = BoxIndexer::UInt(iblock)*MT + threadIdxx;
#endif
             BoxIndexer const& indexer = dp_boxes[ibox];
             if (icell < indexer.numPts()) {
                 IntVect iv = indexer.intVect(icell);
                 for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                     iv[idim] = indexer.lo[idim] + (iv[idim]-indexer.lo[idim])*stride[idim];
                 }
                 auto const d3 = iv.dim3();
                 f(int(ibox), d3.x, d3.y, d3.z);
             }
         });
    AMREX_GPU_ERROR_CHECK();
}

// One launch over the red or black points of all boxes, as in the single-box
// path: thread ic of a row takes i = 2*ic or 2*ic+1, whichever has the color.
template <int MT, FabArrayType MF, typename F>
void
ParallelForRedBlack_doit (MF const& mf, int redblack, F const& f)
{
    const int nboxes = mf.IndexArray().size();
    if (nboxes == 0) { return; }

    auto const& parforinfo = mf.getParForInfoRedBlack();
    const auto nblocks_per_box = parforinfo.getNBlocksPerBox(MT);
    if (nblocks_per_box == 0) { return; }
    AMREX_ASSERT(Long(nblocks_per_box)*Long(nboxes) < Long(std::numeric_limits<int>::max()));
    const int nblocks = nblocks_per_box * nboxes;
    const BoxIndexer* dp_boxes = parforinfo.getBoxes();
    const int* dp_xbounds = parforinfo.getXBounds();
    const Math::FastDivmodU32 fdm_box(nblocks_per_box);

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    amrex::launch_global<MT>
        <<<nblocks, MT, 0, Gpu::gpuStream()>>>
        ([=] AMREX_GPU_DEVICE () noexcept
         {
             std::uint32_t ibox, iblock;
             fdm_box(ibox, iblock, blockIdx.x);
             auto icell = BoxIndexer::UInt(iblock)*MT + threadIdx.x;

#elif defined(AMREX_USE_SYCL)

    amrex::launch<MT>(nblocks, Gpu::gpuStream(),
         [=] AMREX_GPU_DEVICE (sycl::nd_item<1> const& item) noexcept
         {
             auto blockIdxx = std::uint32_t(item.get_group_linear_id());
             int threadIdxx = item.get_local_linear_id();
             std::uint32_t ibox, iblock;
             fdm_box(ibox, iblock, blockIdxx);
             auto icell = BoxIndexer::UInt(iblock)*MT + threadIdxx;
#endif
             // The box coarsened by 2 in x, and the box's x bounds.
             BoxIndexer const& indexer = dp_boxes[ibox];
             if (icell < indexer.numPts()) {
                 auto [ic, j, k] = indexer(icell);
                 const int i = 2*ic + ((j+k+redblack) & 1);
                 if (i >= dp_xbounds[2*ibox] && i <= dp_xbounds[2*ibox+1]) {
                     f(int(ibox), i, j, k);
                 }
             }
         });
    AMREX_GPU_ERROR_CHECK();
}

template <FabArrayType MF, typename F>
void
ParallelForStrided_doit (MF const& mf, IntVect const& stride, IntVect const& offset, F const& f)
{
    const auto& index_array = mf.IndexArray();
    if (index_array.size() == 1) {
        Box const& b = strided_box(mf.box(index_array[0]), stride, offset);
        const IntVect lo = b.smallEnd();
        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            AMREX_D_TERM(i = lo[0] + (i-lo[0])*stride[0];,
                         j = lo[1] + (j-lo[1])*stride[1];,
                         k = lo[2] + (k-lo[2])*stride[2];)
            f(0, i, j, k);
        });
    } else {
#ifdef AMREX_USE_CUDA
        constexpr int MT = 128;
#else
        constexpr int MT = AMREX_GPU_MAX_THREADS;
#endif
        ParallelForStrided_doit<MT>(mf, stride, offset, f);
    }
}

template <FabArrayType MF, typename F>
void
ParallelForRedBlack_doit (MF const& mf, int redblack, F const& f)
{
    const auto& index_array = mf.IndexArray();
    if (index_array.size() == 1) {
        // Checkerboard: i = 2*ic + ((j+k+redblack)&1)
        Box b = mf.box(index_array[0]);
        const int xlo = b.smallEnd(0);
        const int xhi = b.bigEnd(0);
        b.setSmall(0, amrex::coarsen(xlo,2));
        b.setBig  (0, amrex::coarsen(xhi,2));
        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int ic, int j, int k) noexcept
        {
            const int i = 2*ic + ((j+k+redblack) & 1);
            if (i >= xlo && i <= xhi) {
                f(0, i, j, k);
            }
        });
    } else {
#ifdef AMREX_USE_CUDA
        constexpr int MT = 128;
#else
        constexpr int MT = AMREX_GPU_MAX_THREADS;
#endif
        ParallelForRedBlack_doit<MT>(mf, redblack, f);
    }
}

}
/// \endcond

#endif
#endif
