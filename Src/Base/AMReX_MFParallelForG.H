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
void build_par_for_nblocks (char*& a_hp, char*& a_dp, std::pair<int*,int*>& blocks_x, Box*& pboxes,
                            Vector<Box> const& boxes, Vector<Long> const& ncells, int nthreads)
{
    if (!ncells.empty()) {
        const int nboxes = ncells.size();
        const std::size_t nbytes_boxes = amrex::aligned_size(16, (nboxes+1) * sizeof(int));
        const std::size_t nbytes = nbytes_boxes + nboxes*sizeof(Box);
        a_hp = (char*)The_Pinned_Arena()->alloc(nbytes);
        int* hp_blks = (int*)a_hp;
        Box* hp_boxes = (Box*)(a_hp + nbytes_boxes);
        hp_blks[0] = 0;
        Long ntot = 0;
        bool same_size = true;
        for (int i = 0; i < nboxes; ++i) {
            Long nblocks = (ncells[i] + nthreads-1) / nthreads;
            hp_blks[i+1] = hp_blks[i] + static_cast<int>(nblocks);
            ntot += nblocks;
            same_size = same_size && (ncells[i] == ncells[0]);

            new (hp_boxes+i) Box(boxes[i]);
        }
        amrex::ignore_unused(ntot);
        AMREX_ASSERT(static_cast<Long>(hp_blks[nboxes]) == ntot); // no overflow

        a_dp = (char*) The_Arena()->alloc(nbytes);
        Gpu::htod_memcpy_async(a_dp, a_hp, nbytes);

        blocks_x.first = hp_blks;
        blocks_x.second = (same_size) ? nullptr : (int*)a_dp;
        pboxes = (Box*)(a_dp + nbytes_boxes);
    }
}

inline
void destroy_par_for_nblocks (char* hp, char* dp)
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
    auto call_f (F const& f, int b, int i, int j, int k, int n) noexcept
        -> decltype(f(0,0,0,0,0))
    {
        f(b,i,j,k,n);
    }
}

template <int MT, typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& nghost, int ncomp, IntVect const&, bool, F&& f)
{
    const auto& index_array = mf.IndexArray();
    const int nboxes = index_array.size();

    if (nboxes == 0) {
        return;
    } else if (nboxes == 1) {
        Box const& b = amrex::grow(mf.box(index_array[0]), nghost);
        amrex::ParallelFor(b, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            parfor_mf_detail::call_f(f, 0, i, j, k, n);
        });
    } else {
        auto const& parforinfo = mf.getParForInfo(nghost,MT);
        auto par_for_blocks = parforinfo.getBlocks();
        const int nblocks = par_for_blocks.first[nboxes];
        const int block_0_size = par_for_blocks.first[1];
        const int* dp_nblocks = par_for_blocks.second;
        const Box* dp_boxes = parforinfo.getBoxes();

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

        amrex::launch_global<MT>
            <<<nblocks, MT, 0, Gpu::gpuStream()>>>
            ([=] AMREX_GPU_DEVICE () noexcept
             {
                 int ibox, icell;
                 if (dp_nblocks) {
                     ibox = amrex::bisect(dp_nblocks, 0, nboxes, static_cast<int>(blockIdx.x));
                     icell = (blockIdx.x-dp_nblocks[ibox])*MT + threadIdx.x;
                 } else {
                     ibox = blockIdx.x / block_0_size;
                     icell = (blockIdx.x-ibox*block_0_size)*MT + threadIdx.x;
                 }

#elif defined(AMREX_USE_SYCL)

        amrex::launch<MT>(nblocks, Gpu::gpuStream(),
             [=] AMREX_GPU_DEVICE (sycl::nd_item<1> const& item) noexcept
             {
                 int ibox, icell;
                 int blockIdxx = item.get_group_linear_id();
                 int threadIdxx = item.get_local_linear_id();
                 if (dp_nblocks) {
                     ibox = amrex::bisect(dp_nblocks, 0, nboxes, static_cast<int>(blockIdxx));
                     icell = (blockIdxx-dp_nblocks[ibox])*MT + threadIdxx;
                 } else {
                     ibox = blockIdxx / block_0_size;
                     icell = (blockIdxx-ibox*block_0_size)*MT + threadIdxx;
                 }
#endif
                 Box const& b = dp_boxes[ibox];
                 int ncells = b.numPts();
                 if (icell < ncells) {
                     const auto len = amrex::length(b);
                     int k =  icell /   (len.x*len.y);
                     int j = (icell - k*(len.x*len.y)) /   len.x;
                     int i = (icell - k*(len.x*len.y)) - j*len.x;
                     AMREX_D_TERM(i += b.smallEnd(0);,
                                  j += b.smallEnd(1);,
                                  k += b.smallEnd(2);)
                     for (int n = 0; n < ncomp; ++n) {
                         parfor_mf_detail::call_f(f, ibox, i, j, k, n);
                     }
                 }
             });
    }
    AMREX_GPU_ERROR_CHECK();
}

template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& nghost, int ncomp, IntVect const& ts, bool dynamic, F&& f)
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(mf, nghost, ncomp, ts, dynamic, std::forward<F>(f));
}

template <typename MF, typename F>
std::enable_if_t<IsFabArray<MF>::value>
ParallelFor (MF const& mf, IntVect const& nghost, IntVect const& ts, bool dynamic, F&& f)
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(mf, nghost, 1, ts, dynamic, std::forward<F>(f));
}

}

}

#endif
#endif
