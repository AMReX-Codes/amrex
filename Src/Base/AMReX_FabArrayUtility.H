#ifndef AMREX_FABARRAY_UTILITY_H_
#define AMREX_FABARRAY_UTILITY_H_

#include <AMReX_FabArray.H>
#include <AMReX_LayoutData.H>
#include <AMReX_Print.H>
#include <limits>

namespace amrex {

template <class FAB, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceSum (FabArray<FAB> const& fa, int nghost, F f) {
    return ReduceSum(fa, IntVect(nghost), std::move(f));
}

template <class FAB, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceSum (FabArray<FAB> const& fa, IntVect const& nghost, F f)
{
    using value_type = typename FAB::value_type;
    value_type sm = 0;

#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<value_type> ds_sm(sm);
        value_type* d_sm = ds_sm.dataPtr();

        for (MFIter mfi(fa); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB const* fab = fa.fabPtr(mfi);

            const auto ec = amrex::Cuda::ExecutionConfig(bx.numPts()/32L);

            amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(value_type),
                                   Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* block_sum = gsm.dataPtr();
                value_type* sdata = block_sum + 1;

                value_type tsum = 0.0;
                for (auto const tbx : Gpu::Range(bx)) {
                    tsum += f(tbx, *fab);
                }
                sdata[threadIdx.x] = tsum;
                __syncthreads();

                Gpu::blockReduceSum<AMREX_CUDA_MAX_THREADS>(sdata, *block_sum);

                if (threadIdx.x == 0) Cuda::Atomic::Add(d_sm, *block_sum);
            });
        }

        sm = ds_sm.dataValue();
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
           int nghost, F f) {
    return ReduceSum(fa1, fa2, IntVect(nghost), std::move(f));
}

template <class FAB1, class FAB2, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           IntVect const& nghost, F f)
{
    using value_type = typename FAB1::value_type;
    value_type sm = 0;

#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<value_type> ds_sm(sm);
        value_type* d_sm = ds_sm.dataPtr();

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);

            const auto ec = amrex::Cuda::ExecutionConfig(bx.numPts()/32L);

            amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(value_type),
                                   Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* block_sum = gsm.dataPtr();
                value_type* sdata = block_sum + 1;

                value_type tsum = 0.0;
                for (auto const tbx : Gpu::Range(bx)) {
                    tsum += f(tbx, *fab1, *fab2);
                }
                sdata[threadIdx.x] = tsum;
                __syncthreads();

                Gpu::blockReduceSum<AMREX_CUDA_MAX_THREADS>(sdata, *block_sum);

                if (threadIdx.x == 0) Cuda::Atomic::Add(d_sm, *block_sum);
            });
        }

        sm = ds_sm.dataValue();
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
           int nghost, F f) {
  return ReduceSum(fa1, fa2, fa3, nghost, std::move(f));
}

template <class FAB1, class FAB2, class FAB3, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, FabArray<FAB3> const& fa3,
           IntVect const& nghost, F f)
{
    using value_type = typename FAB1::value_type;
    value_type sm = 0;

#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<value_type> ds_sm(sm);
        value_type* d_sm = ds_sm.dataPtr();

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);
            FAB3 const* fab3 = fa3.fabPtr(mfi);

            const auto ec = amrex::Cuda::ExecutionConfig(bx.numPts()/32L);

            amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(value_type),
                                   Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* block_sum = gsm.dataPtr();
                value_type* sdata = block_sum + 1;

                value_type tsum = 0.0;
                for (auto const tbx : Gpu::Range(bx)) {
                    tsum += f(tbx, *fab1, *fab2, *fab3);
                }
                sdata[threadIdx.x] = tsum;
                __syncthreads();

                Gpu::blockReduceSum<AMREX_CUDA_MAX_THREADS>(sdata, *block_sum);

                if (threadIdx.x == 0) Cuda::Atomic::Add(d_sm, *block_sum);
            });
        }

        sm = ds_sm.dataValue();
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
ReduceMin (FabArray<FAB> const& fa, int nghost, F f) {
    return ReduceMin(fa, IntVect(nghost), std::move(f));
}

template <class FAB, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMin (FabArray<FAB> const& fa, IntVect const& nghost, F f) {
    using value_type = typename FAB::value_type;
    constexpr value_type value_max = std::numeric_limits<value_type>::max();
    value_type r = value_max;

#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<value_type> ds_r(r);
        value_type* d_r = ds_r.dataPtr();

        for (MFIter mfi(fa); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB const* fab = fa.fabPtr(mfi);

            const auto ec = amrex::Cuda::ExecutionConfig(bx.numPts()/32L);

            amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(value_type),
                                   Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* block_r = gsm.dataPtr();
                value_type* sdata = block_r + 1;

#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
                value_type tmin = std::numeric_limits<value_type>::max();
#else
                value_type tmin = value_max;
#endif
                for (auto const tbx : Gpu::Range(bx)) {
                    value_type local_tmin = f(tbx, *fab);
                    tmin = amrex::min(tmin, local_tmin);
                }
                sdata[threadIdx.x] = tmin;
                __syncthreads();

                Gpu::blockReduceMin<AMREX_CUDA_MAX_THREADS>(sdata, *block_r);

                if (threadIdx.x == 0) Cuda::Atomic::Min(d_r, *block_r);
            });
        }

        r = ds_r.dataValue();
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
    return ReduceMin(fa1, fa2, IntVect(nghost), std::move(f));
}

template <class FAB1, class FAB2, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMin (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, IntVect const& nghost, F f)
{
    using value_type = typename FAB1::value_type;
    constexpr value_type value_max = std::numeric_limits<value_type>::max();
    value_type r = value_max;

#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<value_type> ds_r(r);
        value_type* d_r = ds_r.dataPtr();

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);

            const auto ec = amrex::Cuda::ExecutionConfig(bx.numPts()/32L);

            amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(value_type),
                                   Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* block_r = gsm.dataPtr();
                value_type* sdata = block_r + 1;

#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
                value_type tmin = std::numeric_limits<value_type>::max();
#else
                value_type tmin = value_max;
#endif
                for (auto const tbx : Gpu::Range(bx)) {
                    value_type local_tmin = f(tbx, *fab1, *fab2);
                    tmin = amrex::min(tmin, local_tmin);
                }
                sdata[threadIdx.x] = tmin;
                __syncthreads();

                Gpu::blockReduceMin<AMREX_CUDA_MAX_THREADS>(sdata, *block_r);

                if (threadIdx.x == 0) Cuda::Atomic::Min(d_r, *block_r);
            });
        }

        r = ds_r.dataValue();
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
    return ReduceMax(fa, IntVect(nghost), std::move(f));
}

template <class FAB, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMax (FabArray<FAB> const& fa, IntVect const& nghost, F f)
{
    using value_type = typename FAB::value_type;
    constexpr value_type value_lowest = std::numeric_limits<value_type>::lowest();
    value_type r = value_lowest;

#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<value_type> ds_r(r);
        value_type* d_r = ds_r.dataPtr();

        for (MFIter mfi(fa); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB const* fab = fa.fabPtr(mfi);

            const auto ec = amrex::Cuda::ExecutionConfig(bx.numPts()/32L);

            amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(value_type),
                                   Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* block_r = gsm.dataPtr();
                value_type* sdata = block_r + 1;

#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
                value_type tmax = std::numeric_limits<value_type>::lowest();
#else
                value_type tmax = value_lowest;
#endif
                for (auto const tbx : Gpu::Range(bx)) {
                    value_type local_tmax = f(tbx, *fab);
                    tmax = amrex::max(tmax, local_tmax);
                }
                sdata[threadIdx.x] = tmax;
                __syncthreads();

                Gpu::blockReduceMax<AMREX_CUDA_MAX_THREADS>(sdata, *block_r);

                if (threadIdx.x == 0) Cuda::Atomic::Max(d_r, *block_r);
            });
        }

        r = ds_r.dataValue();
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
    return ReduceMax(fa1, fa2, IntVect(nghost), std::move(f));
}

template <class FAB1, class FAB2, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMax (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, IntVect const& nghost, F f)
{
    using value_type = typename FAB1::value_type;
    constexpr value_type value_lowest = std::numeric_limits<value_type>::lowest();
    value_type r = value_lowest;

#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<value_type> ds_r(r);
        value_type* d_r = ds_r.dataPtr();

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);

            const auto ec = amrex::Cuda::ExecutionConfig(bx.numPts()/32L);

            amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(value_type),
                                   Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<value_type> gsm;
                value_type* block_r = gsm.dataPtr();
                value_type* sdata = block_r + 1;

#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
                value_type tmax = std::numeric_limits<value_type>::lowest();
#else
                value_type tmax = value_lowest;
#endif
                for (auto const tbx : Gpu::Range(bx)) {
                    value_type local_tmax = f(tbx, *fab1, *fab2);
                    tmax = amrex::max(tmax, local_tmax);
                }
                sdata[threadIdx.x] = tmax;
                __syncthreads();

                Gpu::blockReduceMax<AMREX_CUDA_MAX_THREADS>(sdata, *block_r);

                if (threadIdx.x == 0) Cuda::Atomic::Max(d_r, *block_r);
            });
        }

        r = ds_r.dataValue();
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
    return ReduceLogicalAnd(fa, IntVect(nghost), std::move(f));
}

template <class FAB, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalAnd (FabArray<FAB> const& fa, IntVect const& nghost, F f)
{
    int r = true;

#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<int> ds_r(r);
        int* d_r = ds_r.dataPtr();

        for (MFIter mfi(fa); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB const* fab = fa.fabPtr(mfi);

            const auto ec = amrex::Cuda::ExecutionConfig(bx.numPts()/32L);

            amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(int),
                                   Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<int> gsm;
                int* block_r = gsm.dataPtr();
                int* sdata = block_r + 1;

                int tr = true;
                for (auto const tbx : Gpu::Range(bx)) {
                    tr = tr && f(tbx, *fab);
                }
                sdata[threadIdx.x] = tr;
                __syncthreads();

                Gpu::blockReduceAnd<AMREX_CUDA_MAX_THREADS>(sdata, *block_r);

                if (threadIdx.x == 0) Cuda::Atomic::And(d_r, *block_r);
            });
        }

        r = ds_r.dataValue();
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
    return ReduceLogicalAnd(fa1, fa2, IntVect(nghost), std::move(f));
}

template <class FAB1, class FAB2, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalAnd (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                  IntVect const& nghost, F f)
{
    int r = true;

#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<int> ds_r(r);
        int* d_r = ds_r.dataPtr();

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);

            const auto ec = amrex::Cuda::ExecutionConfig(bx.numPts()/32L);

            amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(int),
                                   Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<int> gsm;
                int* block_r = gsm.dataPtr();
                int* sdata = block_r + 1;

                int tr = true;
                for (auto const tbx : Gpu::Range(bx)) {
                    tr = tr && f(tbx, *fab1, *fab2);
                }
                sdata[threadIdx.x] = tr;
                __syncthreads();

                Gpu::blockReduceAnd<AMREX_CUDA_MAX_THREADS>(sdata, *block_r);

                if (threadIdx.x == 0) Cuda::Atomic::And(d_r, *block_r);
            });
        }

        r = ds_r.dataValue();
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
    return ReduceLogicalOr(fa, IntVect(nghost), std::move(f));
}

template <class FAB, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalOr (FabArray<FAB> const& fa, IntVect const& nghost, F f)
{
    int r = false;

#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<int> ds_r(r);
        int* d_r = ds_r.dataPtr();

        for (MFIter mfi(fa); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB const* fab = fa.fabPtr(mfi);

            const auto ec = amrex::Cuda::ExecutionConfig(bx.numPts()/32L);

            amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(int),
                                   Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<int> gsm;
                int* block_r = gsm.dataPtr();
                int* sdata = block_r + 1;

                int tr = false;
                for (auto const tbx : Gpu::Range(bx)) {
                    tr = tr || f(tbx, *fab);
                }
                sdata[threadIdx.x] = tr;
                __syncthreads();

                Gpu::blockReduceOr<AMREX_CUDA_MAX_THREADS>(sdata, *block_r);

                if (threadIdx.x == 0) Cuda::Atomic::Or(d_r, *block_r);
            });
        }

        r = ds_r.dataValue();
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
    return ReduceLogicalOr(fa1, fa2, IntVect(nghost), std::move(f));
}

template <class FAB1, class FAB2, class F,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalOr (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                 IntVect const& nghost, F f)
{
    int r = false;

#ifdef AMREX_USE_CUDA
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<int> ds_r(r);
        int* d_r = ds_r.dataPtr();

        for (MFIter mfi(fa1); mfi.isValid(); ++mfi)
        {
            const Box& bx = amrex::grow(mfi.validbox(),nghost);
            FAB1 const* fab1 = fa1.fabPtr(mfi);
            FAB2 const* fab2 = fa2.fabPtr(mfi);

            const auto ec = amrex::Cuda::ExecutionConfig(bx.numPts()/32L);

            amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(int),
                                   Gpu::Device::cudaStream()>>>(
            [=] AMREX_GPU_DEVICE () {
                Gpu::SharedMemory<int> gsm;
                int* block_r = gsm.dataPtr();
                int* sdata = block_r + 1;

                int tr = false;
                for (auto const tbx : Gpu::Range(bx)) {
                    tr = tr || f(tbx, *fab1, *fab2);
                }
                sdata[threadIdx.x] = tr;
                __syncthreads();

                Gpu::blockReduceOr<AMREX_CUDA_MAX_THREADS>(sdata, *block_r);

                if (threadIdx.x == 0) Cuda::Atomic::Or(d_r, *block_r);
            });
        }

        r = ds_r.dataValue();
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

template <class FAB, class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
printCell (FabArray<FAB> const& mf, const IntVect& cell, int comp = -1,
           const IntVect& ng = IntVect::TheZeroVector())
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box& bx = amrex::grow(mfi.validbox(), ng);
        if (bx.contains(cell)) {
	    if (comp >= 0) {
                amrex::AllPrint().SetPrecision(17) << " At cell " << cell << " in Box " << bx
                                                   << ": " << mf[mfi](cell, comp) << std::endl;
	    } else {
                std::ostringstream ss;
                ss.precision(17);
                const int ncomp = mf.nComp();
                for (int i = 0; i < ncomp-1; ++i)
		{
                    ss << mf[mfi](cell,i) << ", ";
		}
                ss << mf[mfi](cell,ncomp-1);
                amrex::AllPrint() << " At cell " << cell << " in Box " << bx
                                  << ": " << ss.str() << std::endl;
	    }
        }
    }
}


template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Add (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    Add(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Add (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        if (bx.ok())
        {
            auto const srcFab = src.array(mfi);
            auto       dstFab = dst.array(mfi);
            AMREX_HOST_DEVICE_FOR_4D ( bx, numcomp, i, j, k, n,
            {
                dstFab(i,j,k,n+dstcomp) += srcFab(i,j,k,n+srccomp);
            });
        }
    }
}


template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Copy (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    Copy(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Copy (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        if (bx.ok())
        {
            auto const srcFab = src.array(mfi);
            auto       dstFab = dst.array(mfi);
            AMREX_HOST_DEVICE_FOR_4D ( bx, numcomp, i, j, k, n,
            {
                dstFab(i,j,k,dstcomp+n) = srcFab(i,j,k,srccomp+n);
            });
        }
    }
}


template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Subtract (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    Subtract(dst,src,srccomp,dstcomp,numcomp,nghost);
}

template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Subtract (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        if (bx.ok())
        {
            auto const srcFab = src.array(mfi);
            auto       dstFab = dst.array(mfi);
            AMREX_HOST_DEVICE_FOR_4D ( bx, numcomp, i, j, k, n,
            {
                dstFab(i,j,k,n+dstcomp) -= srcFab(i,j,k,n+srccomp);
            });
        }
    }
}


template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Multiply (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    Multiply(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Multiply (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        if (bx.ok())
        {
            auto const srcFab = src.array(mfi);
            auto       dstFab = dst.array(mfi);
            AMREX_HOST_DEVICE_FOR_4D ( bx, numcomp, i, j, k, n,
            {
                dstFab(i,j,k,n+dstcomp) *= srcFab(i,j,k,n+srccomp);
            });
        }
    }
}


template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Divide (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    Divide(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Divide (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        if (bx.ok())
        {
            auto const srcFab = src.array(mfi);
            auto       dstFab = dst.array(mfi);
            AMREX_HOST_DEVICE_FOR_4D ( bx, numcomp, i, j, k, n,
            {
                dstFab(i,j,k,n+dstcomp) /= srcFab(i,j,k,n+srccomp);
            });
        }
    }
}

template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Abs (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    Abs(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

template <class FAB,
          class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
Abs (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        if (bx.ok())
        {
            auto const srcFab = src.array(mfi);
            auto       dstFab = dst.array(mfi);
            AMREX_HOST_DEVICE_FOR_4D ( bx, numcomp, i, j, k, n,
            {
                dstFab(i,j,k,n+dstcomp) /= srcFab(i,j,k,n+srccomp);
            });
        }
    }
}

template <class FAB, class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
prefetchToHost (FabArray<FAB> const& fa)
{
#ifdef AMREX_USE_CUDA
    for (MFIter mfi(fa); mfi.isValid(); ++mfi) {
        fa.prefetchToHost(mfi);
    }
#endif
}

template <class FAB, class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value> >
void
prefetchToDevice (FabArray<FAB> const& fa)
{
#ifdef AMREX_USE_CUDA
    for (MFIter mfi(fa); mfi.isValid(); ++mfi) {
        fa.prefetchToDevice(mfi);
    }
#endif
}

template <class FAB, class IFAB, class bar = amrex::EnableIf_t<IsBaseFab<FAB>::value
                                                               && IsBaseFab<IFAB>::value> >
void
OverrideSync (FabArray<FAB> & fa, FabArray<IFAB> const& msk, const Periodicity& period)
{
    BL_PROFILE("OverrideSync()");

    if (fa.ixType().cellCentered()) return;
    
    const int ncomp = fa.nComp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fa,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto fab = fa.array(mfi);
        auto const ifab = msk.array(mfi);
        AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            if (!ifab(i,j,k)) fab(i,j,k,n) = 0;
        });
    }
    
    FabArray<FAB> tmpmf(fa.boxArray(), fa.DistributionMap(), ncomp, 0,
                        MFInfo().SetDeviceFab(false), fa.Factory());
    tmpmf.setVal(0);
    tmpmf.ParallelCopy(fa, period, FabArrayBase::ADD);

    amrex::Copy(fa, tmpmf, 0, 0, ncomp, 0);
}

}

#endif
