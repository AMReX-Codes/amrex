#ifndef AMREX_MF_PARALLEL_FOR_G_H_
#define AMREX_MF_PARALLEL_FOR_G_H_
#include <AMReX_Config.H>

#ifdef AMREX_USE_GPU

#include <algorithm>
#include <cmath>
#include <limits>

namespace amrex {
namespace detail {

inline
void build_par_for_boxes (char*& hp, BoxIndexer*& pboxes, Vector<Box> const& boxes)
{
    if (boxes.empty()) { return; }
    const int nboxes = boxes.size();
    const std::size_t nbytes = nboxes*sizeof(BoxIndexer);
    hp = (char*)The_Pinned_Arena()->alloc(nbytes);
    auto* hp_boxes = (BoxIndexer*)hp;
    for (int i = 0; i < nboxes; ++i) {
        new (hp_boxes+i) BoxIndexer(boxes[i]);
    }

    auto dp = (char*) The_Arena()->alloc(nbytes);
    Gpu::htod_memcpy_async(dp, hp, nbytes);
    pboxes = (BoxIndexer*)dp;
}

inline
void destroy_par_for_boxes (char* hp, char* dp)
{
    The_Pinned_Arena()->free(hp);
    The_Arena()->free(dp);
}
}

namespace experimental::detail {

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

template <int MT, typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& nghost, int ncomp, IntVect const&, bool, F const& f)
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

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

        amrex::launch_global<MT>
            <<<nblocks, MT, 0, Gpu::gpuStream()>>>
            ([=] AMREX_GPU_DEVICE () noexcept
             {
                 int ibox = int(blockIdx.x) / nblocks_per_box;
                 auto icell = std::uint64_t(blockIdx.x-ibox*nblocks_per_box)*MT + threadIdx.x;

#elif defined(AMREX_USE_SYCL)

        amrex::launch<MT>(nblocks, Gpu::gpuStream(),
             [=] AMREX_GPU_DEVICE (sycl::nd_item<1> const& item) noexcept
             {
                 int blockIdxx = item.get_group_linear_id();
                 int threadIdxx = item.get_local_linear_id();
                 int ibox = int(blockIdxx) / nblocks_per_box;
                 auto icell = std::uint64_t(blockIdxx-ibox*nblocks_per_box)*MT + threadIdxx;
#endif
                 BoxIndexer const& indexer = dp_boxes[ibox];
                 if (icell < indexer.numPts()) {
                     auto [i, j, k] = indexer(icell);
                     parfor_mf_detail::call_f(f, ibox, i, j, k, ncomp);
                 }
             });
    }
    AMREX_GPU_ERROR_CHECK();
}

template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& nghost, int ncomp, IntVect const& ts, bool dynamic, F&& f)
{
#ifdef AMREX_USE_CUDA
    constexpr int MT = 128;
#else
    constexpr int MT = AMREX_GPU_MAX_THREADS;
#endif
    ParallelFor<MT>(mf, nghost, ncomp, ts, dynamic, std::forward<F>(f));
}

template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& nghost, IntVect const& ts, bool dynamic, F&& f)
{
#ifdef AMREX_USE_CUDA
    constexpr int MT = 128;
#else
    constexpr int MT = AMREX_GPU_MAX_THREADS;
#endif
    ParallelFor<MT>(mf, nghost, 1, ts, dynamic, std::forward<F>(f));
}

}

}

#endif
#endif
