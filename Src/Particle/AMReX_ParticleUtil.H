#ifndef AMREX_PARTICLEUTIL_H_
#define AMREX_PARTICLEUTIL_H_

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_MFIter.H>

#include <limits>

namespace amrex
{

AMREX_GPU_HOST_DEVICE
int getTileIndex (const IntVect& iv, const Box& box, const bool a_do_tiling, 
                  const IntVect& a_tile_size, Box& tbx);


template <class PC, class F>
typename PC::RealType
ReduceRealSum (PC const& pc, int lev, F f)
{
    return ReduceSum(pc, lev, lev, std::move(f));
}

template <class PC, class F>
typename PC::RealType
ReduceSum (PC const& pc, int lev_min, int lev_max, F f)
{
    using value_type = typename PC::RealType;
    value_type sm = 0;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        Gpu::DeviceScalar<value_type> ds_sm(sm);
        value_type* d_sm = ds_sm.dataPtr();

        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            const auto& plevel = pc.GetParticles(lev);
            for(MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
            {
                auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());            
                if(plevel.find(index) == plevel.end()) continue;            
                const auto& tile = plevel.at(index);

                const auto np = tile.numParticles();
                if (np == 0) continue;
            
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();

                constexpr int parts_per_thread = 32;
                const auto ec = amrex::Gpu::ExecutionConfig(bx.numPts()/parts_per_thread);

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
            for(MFIter mfi = pc.MakeMFIter(lev, true); mfi.isValid(); ++mfi)
            {
                auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());            
                if(plevel.find(index) == plevel.end()) continue;            
                const auto& tile = plevel.at(index);

                const auto np = tile.numParticles();
                if (np == 0) continue;
            
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();

                for (int i = 0; i < np; ++i)
                    sm += f(pstruct[i]);
            }
        }
    }

    return sm;
}
    
}

#endif // include guard
