#ifndef AMREX_BASEFAB_UTILITY_H_
#define AMREX_BASEFAB_UTILITY_H_
#include <AMReX_Config.H>

#include <AMReX_BaseFab.H>
#include <AMReX_TableData.H>
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
        constexpr int nthreads_per_block = (STRUCTSIZE <= 8) ? 256 : 128;
        std::uint64_t nblocks_long = (ntotcells+nthreads_per_block-1)/nthreads_per_block;
        AMREX_ASSERT(nblocks_long <= std::uint64_t(std::numeric_limits<int>::max()));
        auto nblocks = int(nblocks_long);
        std::size_t shared_mem_bytes = nthreads_per_block * sizeof(STRUCT);
        T* p = (T*)aos_fab.dataPtr();
#ifdef AMREX_USE_SYCL
        amrex::launch<nthreads_per_block>(nblocks, shared_mem_bytes, Gpu::gpuStream(),
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
        amrex::launch<nthreads_per_block>(nblocks, shared_mem_bytes, Gpu::gpuStream(),
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

//! Transpose 3D array (nx,ny,nz) from row-major (i.e. C order) to
//! column-major (Fortran order). The input's unit stride direction is z,
//! whereas the output's unit stride direction is x.
template <typename T>
void transposeCtoF (T const* pi, T* po, int nx, int ny, int nz)
{
    AMREX_ALWAYS_ASSERT(pi != po);

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    constexpr int tile_dim = 32;
    constexpr int block_rows = 8;
    constexpr int nthreads = tile_dim*block_rows;

    // Each block has tile_dim x block_rows threads. They work on a tile_dim
    // x tile_dim tile.

    dim3 block{unsigned(tile_dim), unsigned(block_rows), 1};
    dim3 grid{unsigned((nx+tile_dim-1)/tile_dim),
              unsigned((nz+tile_dim-1)/tile_dim), unsigned(ny)};

    AMREX_LAUNCH_KERNEL(nthreads, grid, block, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE ()
    {
        __shared__ T tile[tile_dim][tile_dim+1]; // +1 to avoid bank conflicts

        int tx = blockIdx.y * tile_dim + threadIdx.x; // for loading z-direction
        int ty = blockIdx.x * tile_dim + threadIdx.y; // for loading x-direction

        int j = blockIdx.z; // for y-direction

        if (tx < nz) {
            for (int it = 0; it < tile_dim; it += block_rows) {
                int i = it + ty;
                if (i < nx) {
                    //      x               z
                    tile[threadIdx.y+it][threadIdx.x] = pi[tx + (j+i*std::size_t(ny))*nz];
                }
            }
        }

        __syncthreads();

        tx = blockIdx.x * tile_dim + threadIdx.x; // for storing x-direction
        ty = blockIdx.y * tile_dim + threadIdx.y; // for storing z-direction

        if (tx < nx) {
            for (int it = 0; it < tile_dim; it += block_rows) {
                int k = it + ty;
                if (k < nz) {
                    po[tx + (j+k*std::size_t(ny))*nx] = tile[threadIdx.x][threadIdx.y+it];
                }
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();

#else

    Table3D<T const, Order::C> ctab(pi, {0,0,0}, {nx, ny, nz});
    Table3D<T      , Order::F> ftab(po, {0,0,0}, {nx, ny, nz});
    BoxND<3> box(IntVectND<3>(0), IntVectND<3>(nx-1,ny-1,nz-1));
    ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        ftab(i,j,k) = ctab(i,j,k);
    });

#endif

    Gpu::streamSynchronize();
}

//! Transpose 2D array (nx,ny) from row-major (i.e. C order) to column-major
//! (Fortran order). The input's unit stride direction is y, whereas the
//! output's unit stride direction is x.
template <typename T>
void transposeCtoF (T const* pi, T* po, int nx, int ny)
{
    transposeCtoF(pi, po, nx, 1, ny);
}

}

#endif
