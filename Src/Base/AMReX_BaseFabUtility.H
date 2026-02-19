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
//! whereas the output's unit stride direction is x. Note that for GPU
//! builds, the kernel runs on the current GPU stream asynchronously with
//! respect to the host. If synchronization is needed, it's up to the user
//! to call `amrex::Gpu::streamSynchronize()`.
template <typename T>
void transposeCtoF (T const* pi, T* po, int nx, int ny, int nz)
{
    AMREX_ALWAYS_ASSERT(pi != po);

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

    constexpr int tile_dim = 32;
    constexpr int block_rows = 16;
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

        int k = blockIdx.y * tile_dim + threadIdx.x; // for loading z-direction
        int i = blockIdx.x * tile_dim + threadIdx.y; // for loading x-direction

        int j = blockIdx.z; // for y-direction

        if (k < nz) {
            for (int it = 0; it < tile_dim; it += block_rows, i += block_rows) {
                if (i < nx) {
                    //      x               z
                    tile[threadIdx.y+it][threadIdx.x] = pi[k + (j+i*std::size_t(ny))*nz];
                }
            }
        }

        __syncthreads();

        i = blockIdx.x * tile_dim + threadIdx.x; // for storing x-direction
        k = blockIdx.y * tile_dim + threadIdx.y; // for storing z-direction

        if (i < nx) {
            for (int it = 0; it < tile_dim; it += block_rows, k += block_rows) {
                if (k < nz) {
                    po[i + (j+k*std::size_t(ny))*nx] = tile[threadIdx.x][threadIdx.y+it];
                }
            }
        }
    });

#elif defined(AMREX_USE_SYCL)

    constexpr int tile_dim = 32;
    constexpr int block_rows = 8;

    // Each block has tile_dim x block_rows threads. They work on a tile_dim
    // x tile_dim tile.

    sycl::range<3> block{std::size_t(1), std::size_t(block_rows), std::size_t(tile_dim)};
    sycl::range<3> grid{std::size_t(ny), std::size_t((nz+tile_dim-1)/tile_dim),
                                         std::size_t((nx+tile_dim-1)/tile_dim)};
    sycl::range<3> global_size{grid[0]*block[0],
                               grid[1]*block[1],
                               grid[2]*block[2]};

    auto& q = *(Gpu::gpuStream().queue);
    try {
      q.submit([&] (sycl::handler& h)
      {
        auto tile = sycl::local_accessor<T,2>(sycl::range<2>(tile_dim,tile_dim+1),h);

        h.parallel_for(sycl::nd_range<3>(global_size, block),
        [=] (sycl::nd_item<3> item)
        {
            auto group = item.get_group();
            dim3 blockIdx{unsigned(group.get_group_id(2)),
                          unsigned(group.get_group_id(1)),
                          unsigned(group.get_group_id(0))};
            dim3 threadIdx{unsigned(item.get_local_id(2)),
                           unsigned(item.get_local_id(1)),
                           unsigned(item.get_local_id(0))};

            int k = blockIdx.y * tile_dim + threadIdx.x; // for loading z-direction
            int i = blockIdx.x * tile_dim + threadIdx.y; // for loading x-direction

            int j = blockIdx.z; // for y-direction

            if (k < nz) {
                for (int it = 0; it < tile_dim; it += block_rows, i += block_rows) {
                    if (i < nx) {
                        //      x               z
                        tile[threadIdx.y+it][threadIdx.x] = pi[k + (j+i*std::size_t(ny))*nz];
                    }
                }
            }

            item.barrier(sycl::access::fence_space::local_space);

            i = blockIdx.x * tile_dim + threadIdx.x; // for storing x-direction
            k = blockIdx.y * tile_dim + threadIdx.y; // for storing z-direction

            if (i < nx) {
                for (int it = 0; it < tile_dim; it += block_rows, k += block_rows) {
                    if (k < nz) {
                        po[i + (j+k*std::size_t(ny))*nx] = tile[threadIdx.x][threadIdx.y+it];
                    }
                }
            }
        });
      });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("transposeCtoF: ")+ex.what()+"!!!!!");
    }

#else

    constexpr int bx = 32;
    constexpr int bz = 32;

    std::size_t nxy = std::size_t(nx) * ny;
    std::size_t nyz = std::size_t(ny) * nz;

#ifdef AMREX_USE_OMP
#pragma omp parallel for collapse(3)
#endif
    for (int j = 0; j < ny; ++j) {
        for (int k0 = 0; k0 < nz; k0 += bz) {
            for (int i0 = 0; i0 < nx; i0 += bx) {
                int imax = std::min(i0+bx, nx);
                int kmax = std::min(k0+bz, nz);
                auto      * AMREX_RESTRICT pdst = po + j*std::size_t(nx);
                auto const* AMREX_RESTRICT psrc = pi + j*std::size_t(nz);
                for (int i = i0; i < imax; ++i) {
                    AMREX_PRAGMA_SIMD
                    for (int k = k0; k < kmax; ++k) {
                        pdst[i + k*nxy] = psrc[k + i*nyz];
                    }
                }
            }
        }
    }

#endif
}

//! Transpose 2D array (nx,ny) from row-major (i.e. C order) to column-major
//! (Fortran order). The input's unit stride direction is y, whereas the
//! output's unit stride direction is x. Note that for GPU builds, the
//! kernel runs on the current GPU stream asynchronously with respect to the
//! host. If synchronization is needed, it's up to the user to call
//! `amrex::Gpu::streamSynchronize()`.
template <typename T>
void transposeCtoF (T const* pi, T* po, int nx, int ny)
{
    transposeCtoF(pi, po, nx, 1, ny);
}

}

#endif
