#ifndef AMREX_FABARRAY_UTILITY_H_
#define AMREX_FABARRAY_UTILITY_H_

#include <AMReX_FabArray.H>
#include <AMReX_LayoutData.H>
#include <limits>

namespace amrex {

template <class FAB, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceSum (FabArray<FAB> const& fa, int nghost, F f)
{
    using value_type = typename FAB::value_type;
    value_type sm = 0;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        LayoutData<Gpu::GridSize> gs;
        int ntotblocks;
        Gpu::getGridSize(fa,0,gs,ntotblocks);

        Gpu::DeviceVector<value_type> sum_block(ntotblocks);

        for (MFIter mfi(fa); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB const* fab = fa.fabPtr(mfi);

            const auto& grid_size = gs[mfi];
            const int global_bid_begin = grid_size.globalBlockId;
            value_type* pblock = sum_block.dataPtr() + global_bid_begin;

            amrex::launch_global<<<grid_size.numBlocks, grid_size.numThreads,
                grid_size.numThreads*sizeof(value_type),
                Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* sdata = gsm.dataPtr();

                value_type tsum = 0.0;
                for (auto const tbx : Gpu::Range(bx)) {
                    tsum += f(tbx, *fab);
                }
                sdata[threadIdx.x] = tsum;
                __syncthreads();

                Gpu::blockReduceSum<AMREX_CUDA_MAX_THREADS>(sdata, *(pblock+blockIdx.x));
            });
        }

        sm = thrust::reduce(sum_block.begin(), sum_block.end(), 0.0, thrust::plus<value_type>());
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(fa,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            sm += f(bx, fa[mfi]);
        }
    }

    return sm;
}

template <class FAB1, class FAB2, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           int nghost, F f)
{
    using value_type = typename FAB1::value_type;
    value_type sm = 0;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        LayoutData<Gpu::GridSize> gs;
        int ntotblocks;
        Gpu::getGridSize(fa1,0,gs,ntotblocks);

        Gpu::DeviceVector<value_type> sum_block(ntotblocks);

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);

            const auto& grid_size = gs[mfi];
            const int global_bid_begin = grid_size.globalBlockId;
            value_type* pblock = sum_block.dataPtr() + global_bid_begin;

            amrex::launch_global<<<grid_size.numBlocks, grid_size.numThreads,
                grid_size.numThreads*sizeof(value_type),
                Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* sdata = gsm.dataPtr();

                value_type tsum = 0.0;
                for (auto const tbx : Gpu::Range(bx)) {
                    tsum += f(tbx, *fab1, *fab2);
                }
                sdata[threadIdx.x] = tsum;
                __syncthreads();

                Gpu::blockReduceSum<AMREX_CUDA_MAX_THREADS>(sdata, *(pblock+blockIdx.x));
            });
        }

        sm = thrust::reduce(sum_block.begin(), sum_block.end(), 0.0, thrust::plus<value_type>());
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            sm += f(bx, fa1[mfi], fa2[mfi]);
        }
    }

    return sm;
}

template <class FAB1, class FAB2, class FAB3, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, FabArray<FAB3> const& fa3,
           int nghost, F f)
{
    using value_type = typename FAB1::value_type;
    value_type sm = 0;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        LayoutData<Gpu::GridSize> gs;
        int ntotblocks;
        Gpu::getGridSize(fa1,0,gs,ntotblocks);

        Gpu::DeviceVector<value_type> sum_block(ntotblocks);

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);
            FAB3 const* fab3 = fa3.fabPtr(mfi);

            const auto& grid_size = gs[mfi];
            const int global_bid_begin = grid_size.globalBlockId;
            value_type* pblock = sum_block.dataPtr() + global_bid_begin;

            amrex::launch_global<<<grid_size.numBlocks, grid_size.numThreads,
                grid_size.numThreads*sizeof(value_type),
                Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* sdata = gsm.dataPtr();

                value_type tsum = 0.0;
                for (auto const tbx : Gpu::Range(bx)) {
                    tsum += f(tbx, *fab1, *fab2, *fab3);
                }
                sdata[threadIdx.x] = tsum;
                __syncthreads();

                Gpu::blockReduceSum<AMREX_CUDA_MAX_THREADS>(sdata, *(pblock+blockIdx.x));
            });
        }

        sm = thrust::reduce(sum_block.begin(), sum_block.end(), 0.0, thrust::plus<value_type>());
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            sm += f(bx, fa1[mfi], fa2[mfi], fa3[mfi]);
        }
    }

    return sm;
}

template <class FAB, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMin (FabArray<FAB> const& fa, int nghost, F f)
{
    using value_type = typename FAB::value_type;
    value_type r = std::numeric_limits<value_type>::max();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        LayoutData<Gpu::GridSize> gs;
        int ntotblocks;
        Gpu::getGridSize(fa,0,gs,ntotblocks);

        Gpu::DeviceVector<value_type> tmin_block(ntotblocks);

        for (MFIter mfi(fa); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB const* fab = fa.fabPtr(mfi);

            const auto& grid_size = gs[mfi];
            const int global_bid_begin = grid_size.globalBlockId;
            value_type* pblock = tmin_block.dataPtr() + global_bid_begin;

            amrex::launch_global<<<grid_size.numBlocks, grid_size.numThreads,
                grid_size.numThreads*sizeof(value_type), Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* sdata = gsm.dataPtr();

                value_type tmin = std::numeric_limits<value_type>::max();
                for (auto const tbx : Gpu::Range(bx)) {
                    value_type local_tmin = f(tbx, *fab);
                    tmin = amrex::min(tmin, local_tmin);
                }
                sdata[threadIdx.x] = tmin;
                __syncthreads();

                Gpu::blockReduceMin<AMREX_CUDA_MAX_THREADS>(sdata, *(pblock+blockIdx.x));
            });
        }

        r = *(thrust::min_element(tmin_block.begin(), tmin_block.end()));
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(min:r)
#endif
        for (MFIter mfi(fa,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            r = std::min(r, f(bx, fa[mfi]));
        }
    }

    return r;
}

template <class FAB1, class FAB2, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMin (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, int nghost, F f)
{
    using value_type = typename FAB1::value_type;
    value_type r = std::numeric_limits<value_type>::max();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        LayoutData<Gpu::GridSize> gs;
        int ntotblocks;
        Gpu::getGridSize(fa1,0,gs,ntotblocks);

        Gpu::DeviceVector<value_type> tmin_block(ntotblocks);

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);

            const auto& grid_size = gs[mfi];
            const int global_bid_begin = grid_size.globalBlockId;
            value_type* pblock = tmin_block.dataPtr() + global_bid_begin;

            amrex::launch_global<<<grid_size.numBlocks, grid_size.numThreads,
                grid_size.numThreads*sizeof(value_type), Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* sdata = gsm.dataPtr();

                value_type tmin = std::numeric_limits<value_type>::max();
                for (auto const tbx : Gpu::Range(bx)) {
                    value_type local_tmin = f(tbx, *fab1, *fab2);
                    tmin = amrex::min(tmin, local_tmin);
                }
                sdata[threadIdx.x] = tmin;
                __syncthreads();

                Gpu::blockReduceMin<AMREX_CUDA_MAX_THREADS>(sdata, *(pblock+blockIdx.x));
            });
        }

        r = *(thrust::min_element(tmin_block.begin(), tmin_block.end()));
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(min:r)
#endif
        for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            r = std::min(r, f(bx, fa1[mfi], fa2[mfi]));
        }
    }

    return r;
}

template <class FAB, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMax (FabArray<FAB> const& fa, int nghost, F f)
{
    using value_type = typename FAB::value_type;
    value_type r = std::numeric_limits<value_type>::lowest();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        LayoutData<Gpu::GridSize> gs;
        int ntotblocks;
        Gpu::getGridSize(fa,0,gs,ntotblocks);

        Gpu::DeviceVector<value_type> tmax_block(ntotblocks);

        for (MFIter mfi(fa); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB const* fab = fa.fabPtr(mfi);

            const auto& grid_size = gs[mfi];
            const int global_bid_begin = grid_size.globalBlockId;
            value_type* pblock = tmax_block.dataPtr() + global_bid_begin;

            amrex::launch_global<<<grid_size.numBlocks, grid_size.numThreads,
                grid_size.numThreads*sizeof(value_type), Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* sdata = gsm.dataPtr();

                value_type tmax = std::numeric_limits<value_type>::lowest();
                for (auto const tbx : Gpu::Range(bx)) {
                    value_type local_tmax = f(tbx, *fab);
                    tmax = amrex::max(tmax, local_tmax);
                }
                sdata[threadIdx.x] = tmax;
                __syncthreads();

                Gpu::blockReduceMax<AMREX_CUDA_MAX_THREADS>(sdata, *(pblock+blockIdx.x));
            });
        }

        r = *(thrust::max_element(tmax_block.begin(), tmax_block.end()));
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(max:r)
#endif
        for (MFIter mfi(fa,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            r = std::max(r, f(bx, fa[mfi]));
        }
    }

    return r;
}

template <class FAB1, class FAB2, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMax (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, int nghost, F f)
{
    using value_type = typename FAB1::value_type;
    value_type r = std::numeric_limits<value_type>::lowest();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        LayoutData<Gpu::GridSize> gs;
        int ntotblocks;
        Gpu::getGridSize(fa1,0,gs,ntotblocks);

        Gpu::DeviceVector<value_type> tmax_block(ntotblocks);

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);

            const auto& grid_size = gs[mfi];
            const int global_bid_begin = grid_size.globalBlockId;
            value_type* pblock = tmax_block.dataPtr() + global_bid_begin;

            amrex::launch_global<<<grid_size.numBlocks, grid_size.numThreads,
                grid_size.numThreads*sizeof(value_type), Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* sdata = gsm.dataPtr();

                value_type tmax = std::numeric_limits<value_type>::lowest();
                for (auto const tbx : Gpu::Range(bx)) {
                    value_type local_tmax = f(tbx, *fab1, *fab2);
                    tmax = amrex::max(tmax, local_tmax);
                }
                sdata[threadIdx.x] = tmax;
                __syncthreads();

                Gpu::blockReduceMax<AMREX_CUDA_MAX_THREADS>(sdata, *(pblock+blockIdx.x));
            });
        }

        r = *(thrust::max_element(tmax_block.begin(), tmax_block.end()));
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(max:r)
#endif
        for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            r = std::max(r, f(bx, fa1[mfi], fa2[mfi]));
        }
    }

    return r;
}

template <class FAB, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalAnd (FabArray<FAB> const& fa, int nghost, F f)
{
    bool r = true;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        LayoutData<Gpu::GridSize> gs;
        int ntotblocks;
        Gpu::getGridSize(fa,0,gs,ntotblocks);

        Gpu::DeviceVector<unsigned char> tr_block(ntotblocks);

        for (MFIter mfi(fa); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB const* fab = fa.fabPtr(mfi);

            const auto& grid_size = gs[mfi];
            const int global_bid_begin = grid_size.globalBlockId;
            unsigned char* pblock = tr_block.dataPtr() + global_bid_begin;

            amrex::launch_global<<<grid_size.numBlocks, grid_size.numThreads,
                grid_size.numThreads*sizeof(unsigned char), Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<unsigned char> gsm;
                unsigned char* sdata = gsm.dataPtr();

                bool tr = true;
                for (auto const tbx : Gpu::Range(bx)) {
                    tr = tr && f(tbx, *fab);
                }
                sdata[threadIdx.x] = tr;
                __syncthreads();

                Gpu::blockReduceAnd<AMREX_CUDA_MAX_THREADS>(sdata, *(pblock+blockIdx.x));
            });
        }

        r = thrust::reduce(tr_block.begin(), tr_block.end(), true, thrust::logical_and<unsigned char>());
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(&&:r)
#endif
        for (MFIter mfi(fa,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            r = r && f(bx, fa[mfi]);
        }
    }

    return r;
}

template <class FAB1, class FAB2, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalAnd (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                  int nghost, F f)
{
    bool r = true;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        LayoutData<Gpu::GridSize> gs;
        int ntotblocks;
        Gpu::getGridSize(fa1,0,gs,ntotblocks);

        Gpu::DeviceVector<unsigned char> tr_block(ntotblocks);

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);

            const auto& grid_size = gs[mfi];
            const int global_bid_begin = grid_size.globalBlockId;
            unsigned char* pblock = tr_block.dataPtr() + global_bid_begin;

            amrex::launch_global<<<grid_size.numBlocks, grid_size.numThreads,
                grid_size.numThreads*sizeof(unsigned char), Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<unsigned char> gsm;
                unsigned char* sdata = gsm.dataPtr();

                bool tr = true;
                for (auto const tbx : Gpu::Range(bx)) {
                    tr = tr && f(tbx, *fab1, *fab2);
                }
                sdata[threadIdx.x] = tr;
                __syncthreads();

                Gpu::blockReduceAnd<AMREX_CUDA_MAX_THREADS>(sdata, *(pblock+blockIdx.x));
            });
        }

        r = thrust::reduce(tr_block.begin(), tr_block.end(), true, thrust::logical_and<unsigned char>());
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(&&:r)
#endif
        for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            r = r && f(bx, fa1[mfi], fa2[mfi]);
        }
    }

    return r;
}

template <class FAB, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalOr (FabArray<FAB> const& fa, int nghost, F f)
{
    bool r = false;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        LayoutData<Gpu::GridSize> gs;
        int ntotblocks;
        Gpu::getGridSize(fa,0,gs,ntotblocks);

        Gpu::DeviceVector<unsigned char> tr_block(ntotblocks);

        for (MFIter mfi(fa); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB const* fab = fa.fabPtr(mfi);

            const auto& grid_size = gs[mfi];
            const int global_bid_begin = grid_size.globalBlockId;
            unsigned char* pblock = tr_block.dataPtr() + global_bid_begin;

            amrex::launch_global<<<grid_size.numBlocks, grid_size.numThreads,
                grid_size.numThreads*sizeof(unsigned char), Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<unsigned char> gsm;
                unsigned char* sdata = gsm.dataPtr();

                bool tr = false;
                for (auto const tbx : Gpu::Range(bx)) {
                    tr = tr || f(tbx, *fab);
                }
                sdata[threadIdx.x] = tr;
                __syncthreads();

                Gpu::blockReduceOr<AMREX_CUDA_MAX_THREADS>(sdata, *(pblock+blockIdx.x));
            });
        }

        r = thrust::reduce(tr_block.begin(), tr_block.end(), false, thrust::logical_or<unsigned char>());
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(||:r)
#endif
        for (MFIter mfi(fa,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            r = r || f(bx, fa[mfi]);
        }
    }

    return r;
}

template <class FAB1, class FAB2, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalOr (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                 int nghost, F f)
{
    bool r = false;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        LayoutData<Gpu::GridSize> gs;
        int ntotblocks;
        Gpu::getGridSize(fa1,0,gs,ntotblocks);

        Gpu::DeviceVector<unsigned char> tr_block(ntotblocks);

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);

            const auto& grid_size = gs[mfi];
            const int global_bid_begin = grid_size.globalBlockId;
            unsigned char* pblock = tr_block.dataPtr() + global_bid_begin;

            amrex::launch_global<<<grid_size.numBlocks, grid_size.numThreads,
                grid_size.numThreads*sizeof(unsigned char), Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<unsigned char> gsm;
                unsigned char* sdata = gsm.dataPtr();

                bool tr = false;
                for (auto const tbx : Gpu::Range(bx)) {
                    tr = tr || f(tbx, *fab1, *fab2);
                }
                sdata[threadIdx.x] = tr;
                __syncthreads();

                Gpu::blockReduceOr<AMREX_CUDA_MAX_THREADS>(sdata, *(pblock+blockIdx.x));
            });
        }

        r = thrust::reduce(tr_block.begin(), tr_block.end(), false, thrust::logical_or<unsigned char>());
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(||:r)
#endif
        for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            r = r || f(bx, fa1[mfi], fa2[mfi]);
        }
    }

    return r;
}

}

#endif
