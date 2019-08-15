#ifndef AMREX_PARTICLEREDUCE_H_
#define AMREX_PARTICLEREDUCE_H_

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_GpuUtility.H>

#include <limits>

namespace amrex
{

template <class PC, class F>
auto
ReduceSum (PC const& pc, F f) -> decltype(f(typename PC::ParticleType()))
{
    return ReduceSum(pc, 0, pc.finestLevel(), std::move(f));
}
    
template <class PC, class F>
auto
ReduceSum (PC const& pc, int lev, F f) -> decltype(f(typename PC::ParticleType()))
{
    return ReduceSum(pc, lev, lev, std::move(f));
}

template <class PC, class F>
auto
ReduceSum (PC const& pc, int lev_min, int lev_max, F f) -> decltype(f(typename PC::ParticleType()))
{
    using value_type = decltype(f(typename PC::ParticleType()));
    using ParIter = typename PC::ParConstIterType;
    value_type sm = 0;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<value_type> ds_sm(sm);
        value_type* d_sm = ds_sm.dataPtr();

        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();        

                constexpr int parts_per_thread = 32;
                const auto ec = amrex::Gpu::ExecutionConfig(np/parts_per_thread);

                amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(value_type),
                                   Gpu::gpuStream()>>>(
                [=] AMREX_GPU_DEVICE () {
                    Gpu::SharedMemory<value_type> gsm;
                    value_type* block_sum = gsm.dataPtr();
                    value_type* sdata = block_sum + 1;

                    value_type tsum = 0.0;
                    for (auto const i : Gpu::Range(np)) {
                        tsum += f(pstruct[i]);
                    }
                    sdata[threadIdx.x] = tsum;
                    __syncthreads();

                    Gpu::blockReduceSum<AMREX_GPU_MAX_THREADS,Gpu::Device::warp_size>(sdata, *block_sum);

                    if (threadIdx.x == 0) Gpu::Atomic::Add(d_sm, *block_sum);
                });
            }
        }
        sm = ds_sm.dataValue();
    }
    else
#endif
    {
        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();
                for (int i = 0; i < np; ++i)
                    sm += f(pstruct[i]);
            }
        }
    }

    return sm;
}

template <class PC, class F>
auto
ReduceMax (PC const& pc, F f) -> decltype(f(typename PC::ParticleType()))
{
    return ReduceMax(pc, 0, pc.finestLevel(), std::move(f));
}
    
template <class PC, class F>
auto
ReduceMax (PC const& pc, int lev, F f) -> decltype(f(typename PC::ParticleType()))
{
    return ReduceMax(pc, lev, lev, std::move(f));
}

template <class PC, class F>
auto
ReduceMax (PC const& pc, int lev_min, int lev_max, F f) -> decltype(f(typename PC::ParticleType()))
{
    using value_type = decltype(f(typename PC::ParticleType()));
    using ParIter = typename PC::ParConstIterType;
    constexpr value_type value_lowest = std::numeric_limits<value_type>::lowest();
    value_type r = value_lowest;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<value_type> ds_r(r);
        value_type* d_r = ds_r.dataPtr();

        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();        

                constexpr int parts_per_thread = 32;
                const auto ec = amrex::Gpu::ExecutionConfig(np/parts_per_thread);

                amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(value_type),
                                   Gpu::gpuStream()>>>(
                [=] AMREX_GPU_DEVICE () {
                    Gpu::SharedMemory<value_type> gr;
                    value_type* block_r = gr.dataPtr();
                    value_type* sdata = block_r + 1;

                    value_type tmax = value_lowest;
                    for (auto const i : Gpu::Range(np)) {
                        value_type local_tmax = f(pstruct[i]);
                        tmax = amrex::max(tmax, local_tmax);
                    }
                    sdata[threadIdx.x] = tmax;
                    __syncthreads();

                    Gpu::blockReduceMax<AMREX_GPU_MAX_THREADS,Gpu::Device::warp_size>(sdata, *block_r);

                    if (threadIdx.x == 0) Gpu::Atomic::Max(d_r, *block_r);
                });
            }
        }
        r = ds_r.dataValue();
    }
    else
#endif
    {
        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(max:r)
#endif
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();        

                for (int i = 0; i < np; ++i)
                    r = std::max(r, f(pstruct[i]));
            }
        }
    }

    return r;
}

template <class PC, class F>
auto
ReduceMin (PC const& pc, F f) -> decltype(f(typename PC::ParticleType()))
{
    return ReduceMin(pc, 0, pc.finestLevel(), std::move(f));
}

template <class PC, class F>
auto
ReduceMin (PC const& pc, int lev, F f) -> decltype(f(typename PC::ParticleType()))
{
    return ReduceMin(pc, lev, lev, std::move(f));
}

template <class PC, class F>
auto
ReduceMin (PC const& pc, int lev_min, int lev_max, F f) -> decltype(f(typename PC::ParticleType()))
{
    using value_type = decltype(f(typename PC::ParticleType()));
    using ParIter = typename PC::ParConstIterType;
    constexpr value_type value_max = std::numeric_limits<value_type>::max();
    value_type r = value_max;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<value_type> ds_r(r);
        value_type* d_r = ds_r.dataPtr();

        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();        

                constexpr int parts_per_thread = 32;
                const auto ec = amrex::Gpu::ExecutionConfig(np/parts_per_thread);

                amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(value_type),
                                   Gpu::gpuStream()>>>(
                [=] AMREX_GPU_DEVICE () {
                    Gpu::SharedMemory<value_type> gr;
                    value_type* block_r = gr.dataPtr();
                    value_type* sdata = block_r + 1;

                    value_type tmin = value_max;
                    for (auto const i : Gpu::Range(np)) {
                        value_type local_tmin = f(pstruct[i]);
                        tmin = amrex::min(tmin, local_tmin);
                    }
                    sdata[threadIdx.x] = tmin;
                    __syncthreads();

                    Gpu::blockReduceMin<AMREX_GPU_MAX_THREADS,Gpu::Device::warp_size>(sdata, *block_r);

                    if (threadIdx.x == 0) Gpu::Atomic::Min(d_r, *block_r);
                });
            }
        }
        r = ds_r.dataValue();
    }
    else
#endif
    {
        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(min:r)
#endif
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();        

                for (int i = 0; i < np; ++i)
                    r = std::min(r, f(pstruct[i]));
            }
        }
    }

    return r;
}

template <class PC, class F>
bool
ReduceLogicalAnd (PC const& pc, F f)
{
    return ReduceLogicalAnd(pc, 0, pc.finestLevel(), std::move(f));
}

template <class PC, class F>
bool
ReduceLogicalAnd (PC const& pc, int lev, F f)
{
    return ReduceLogicalAnd(pc, lev, lev, std::move(f));
}

template <class PC, class F>
bool
ReduceLogicalAnd (PC const& pc, int lev_min, int lev_max, F f)
{
    using ParIter = typename PC::ParConstIterType;
    int r = true;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<int> ds_r(r);
        int* d_r = ds_r.dataPtr();

        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();        

                constexpr int parts_per_thread = 32;
                const auto ec = amrex::Gpu::ExecutionConfig(np/parts_per_thread);

                amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(int),
                                   Gpu::gpuStream()>>>(
                [=] AMREX_GPU_DEVICE () {
                    Gpu::SharedMemory<int> gr;
                    int* block_r = gr.dataPtr();
                    int* sdata = block_r + 1;

                    int tr = true;
                    for (auto const i : Gpu::Range(np)) {
                        tr = tr && f(pstruct[i]);
                    }
                    sdata[threadIdx.x] = tr;
                    __syncthreads();

                    Gpu::blockReduceAnd<AMREX_GPU_MAX_THREADS,Gpu::Device::warp_size>(sdata, *block_r);

                    if (threadIdx.x == 0) Gpu::Atomic::LogicalAnd(d_r, *block_r);
                });
            }
        }
        r = ds_r.dataValue();
    }
    else
#endif
    {
        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(min:r)
#endif
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();        

                for (int i = 0; i < np; ++i)
                    r = r && f(pstruct[i]);
            }
        }
    }

    return r;
}

template <class PC, class F>
bool
ReduceLogicalOr (PC const& pc, F f)
{
    return ReduceLogicalOr(pc, 0, pc.finestLevel(), std::move(f));
}

template <class PC, class F>
bool
ReduceLogicalOr (PC const& pc, int lev, F f)
{
    return ReduceLogicalOr(pc, lev, lev, std::move(f));
}

template <class PC, class F>
bool
ReduceLogicalOr (PC const& pc, int lev_min, int lev_max, F f)
{
    using ParIter = typename PC::ParConstIterType;
    int r = false;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<int> ds_r(r);
        int* d_r = ds_r.dataPtr();

        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();        

                constexpr int parts_per_thread = 32;
                const auto ec = amrex::Gpu::ExecutionConfig(np/parts_per_thread);

                amrex::launch_global<<<ec.numBlocks, ec.numThreads, (ec.numThreads.x+1)*sizeof(int),
                                   Gpu::gpuStream()>>>(
                [=] AMREX_GPU_DEVICE () {
                    Gpu::SharedMemory<int> gr;
                    int* block_r = gr.dataPtr();
                    int* sdata = block_r + 1;

                    int tr = false;
                    for (auto const i : Gpu::Range(np)) {
                        tr = tr || f(pstruct[i]);
                    }
                    sdata[threadIdx.x] = tr;
                    __syncthreads();

                    Gpu::blockReduceOr<AMREX_GPU_MAX_THREADS,Gpu::Device::warp_size>(sdata, *block_r);

                    if (threadIdx.x == 0) Gpu::Atomic::LogicalOr(d_r, *block_r);
                });
            }
        }
        r = ds_r.dataValue();
    }
    else
#endif
    {
        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
#ifdef _OPENMP
#pragma omp parallel if (!system::regtest_reduction) reduction(min:r)
#endif
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();        

                for (int i = 0; i < np; ++i)
                    r = r || f(pstruct[i]);
            }
        }
    }

    return r;
}    
}
#endif
