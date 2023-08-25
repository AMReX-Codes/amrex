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
          typename std::enable_if<(sizeof(STRUCT)<=36*8) &&
                                  AMREX_IS_TRIVIALLY_COPYABLE(STRUCT) &&
                                  std::is_trivially_destructible<STRUCT>::value,
                                  int>::type FOO = 0>
void fill (BaseFab<STRUCT>& aos_fab, F && f)
{
    Box const& box = aos_fab.box();
    auto const& aos = aos_fab.array();
    using T = typename STRUCT::value_type;
    constexpr int STRUCTSIZE = sizeof(STRUCT)/sizeof(T);
    static_assert(sizeof(STRUCT) == sizeof(T)*STRUCTSIZE,
                  "amrex::fill: sizeof(STRUCT) != sizeof(T)*STRUCTSIZE");
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        const auto lo  = amrex::lbound(box);
        const auto len = amrex::length(box);
        const auto lenxy = len.x*len.y;
        const auto lenx = len.x;
        int ntotcells = box.numPts();
        int nthreads_per_block = (STRUCTSIZE <= 8) ? 256 : 128;
        int nblocks = (ntotcells+nthreads_per_block-1)/nthreads_per_block;
        std::size_t shared_mem_bytes = nthreads_per_block * sizeof(STRUCT);
        T* p = (T*)aos_fab.dataPtr();
#ifdef AMREX_USE_SYCL
        amrex::launch(nblocks, nthreads_per_block, shared_mem_bytes, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& handler) noexcept
        {
            int icell = handler.globalIdx();
            unsigned int blockDimx = handler.blockDim();
            unsigned int threadIdxx = handler.threadIdx();
            unsigned int blockIdxx = handler.blockIdx();
            auto const shared = (T*)handler.sharedMemory();
            if (icell < ntotcells) {
                auto ga = new(shared+threadIdxx*STRUCTSIZE) STRUCT;
                int k =  icell /   lenxy;
                int j = (icell - k*lenxy) /   lenx;
                int i = (icell - k*lenxy) - j*lenx;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                f(*ga, i, j, k);
            }
            handler.sharedBarrier();
            for (unsigned int m = threadIdxx,
                     mend = amrex::min<unsigned int>(blockDimx, ntotcells-blockDimx*blockIdxx) * STRUCTSIZE;
                 m < mend; m += blockDimx) {
                p[blockDimx*blockIdxx*STRUCTSIZE+m] = shared[m];
            }
        });
#else
        amrex::launch(nblocks, nthreads_per_block, shared_mem_bytes, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept
        {
            int icell = blockDim.x*blockIdx.x+threadIdx.x;
            Gpu::SharedMemory<T> gsm;
            T* const shared = gsm.dataPtr();
            if (icell < ntotcells) {
                auto ga = new(shared+threadIdx.x*STRUCTSIZE) STRUCT;
                int k =  icell /   lenxy;
                int j = (icell - k*lenxy) /   lenx;
                int i = (icell - k*lenxy) - j*lenx;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                f(*ga, i, j, k);
            }
            __syncthreads();
            for (unsigned int m = threadIdx.x,
                     mend = amrex::min<unsigned int>(blockDim.x, ntotcells-blockDim.x*blockIdx.x) * STRUCTSIZE;
                 m < mend; m += blockDim.x) {
                p[blockDim.x*blockIdx.x*STRUCTSIZE+m] = shared[m];
            }
        });
#endif
    } else
#endif
    {
        amrex::LoopOnCpu(box, [=] (int i, int j, int k) noexcept
        {
            f(aos(i,j,k), i, j, k);
        });
    }
}

}

#endif
