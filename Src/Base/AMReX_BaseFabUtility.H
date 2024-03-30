#ifndef AMREX_BASEFAB_UTILITY_H_
#define AMREX_BASEFAB_UTILITY_H_
#include <AMReX_Config.H>

#include <AMReX_BaseFab.H>
#include <AMReX_TypeTraits.H>

namespace amrex {

template <class Tto, class Tfrom>
AMREX_GPU_HOST_DEVICE
void
cast (BaseFab<Tto>& tofab, BaseFab<Tfrom> const& fromfab,
      Box const& bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    auto const& tdata = tofab.array();
    auto const& fdata = fromfab.const_array();
    amrex::LoopConcurrent(bx, ncomp.n, [=] (int i, int j, int k, int n) noexcept
    {
        tdata(i,j,k,n+dcomp.i) = static_cast<Tto>(fdata(i,j,k,n+scomp.i));
    });
}

template <typename STRUCT, typename F,
          std::enable_if_t<(sizeof(STRUCT)<=36*8) &&
                                  AMREX_IS_TRIVIALLY_COPYABLE(STRUCT) &&
                                  std::is_trivially_destructible_v<STRUCT>,
                                  int>FOO = 0>
void fill (BaseFab<STRUCT>& aos_fab, F const& f)
{
    Box const& box = aos_fab.box();
    auto const& aos = aos_fab.array();
    using T = typename STRUCT::value_type;
    constexpr int STRUCTSIZE = sizeof(STRUCT)/sizeof(T);
    static_assert(sizeof(STRUCT) == sizeof(T)*STRUCTSIZE,
                  "amrex::fill: sizeof(STRUCT) != sizeof(T)*STRUCTSIZE");
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        BoxIndexer indexer(box);
        const auto ntotcells = std::uint64_t(box.numPts());
        int nthreads_per_block = (STRUCTSIZE <= 8) ? 256 : 128;
        std::uint64_t nblocks_long = (ntotcells+nthreads_per_block-1)/nthreads_per_block;
        AMREX_ASSERT(nblocks_long <= std::uint64_t(std::numeric_limits<int>::max()));
        auto nblocks = int(nblocks_long);
        std::size_t shared_mem_bytes = nthreads_per_block * sizeof(STRUCT);
        T* p = (T*)aos_fab.dataPtr();
#ifdef AMREX_USE_SYCL
        amrex::launch(nblocks, nthreads_per_block, shared_mem_bytes, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& handler) noexcept
        {
            auto const icell = std::uint64_t(handler.globalIdx());
            std::uint64_t const blockDimx = handler.blockDim();
            std::uint64_t const threadIdxx = handler.threadIdx();
            std::uint64_t const blockIdxx = handler.blockIdx();
            auto const shared = (T*)handler.sharedMemory();
            if (icell < indexer.numPts()) {
                auto ga = new(shared+threadIdxx*STRUCTSIZE) STRUCT;
                auto [i, j, k] = indexer(icell);
                f(*ga, i, j, k);
            }
            handler.sharedBarrier();
            for (std::uint64_t m = threadIdxx,
                     mend = amrex::min<std::uint64_t>(blockDimx, indexer.numPts()-blockDimx*blockIdxx) * STRUCTSIZE;
                 m < mend; m += blockDimx) {
                p[blockDimx*blockIdxx*STRUCTSIZE+m] = shared[m];
            }
        });
#else
        amrex::launch(nblocks, nthreads_per_block, shared_mem_bytes, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept
        {
            std::uint64_t const icell = std::uint64_t(blockDim.x)*blockIdx.x+threadIdx.x;
            Gpu::SharedMemory<T> gsm;
            T* const shared = gsm.dataPtr();
            if (icell < indexer.numPts()) {
                auto ga = new(shared+std::uint64_t(threadIdx.x)*STRUCTSIZE) STRUCT;
                auto [i, j, k] = indexer(icell);
                f(*ga, i, j, k);
            }
            __syncthreads();
            for (std::uint64_t m = threadIdx.x,
                     mend = amrex::min<std::uint64_t>(blockDim.x, indexer.numPts()-std::uint64_t(blockDim.x)*blockIdx.x) * STRUCTSIZE;
                 m < mend; m += blockDim.x) {
                p[std::uint64_t(blockDim.x)*blockIdx.x*STRUCTSIZE+m] = shared[m];
            }
        });
#endif
    } else
#endif
    {
        amrex::LoopOnCpu(box, [&] (int i, int j, int k) noexcept
        {
            f(aos(i,j,k), i, j, k);
        });
    }
}

}

#endif
