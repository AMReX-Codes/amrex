#ifndef AMREX_PARTICLEMESH_H_
#define AMREX_PARTICLEMESH_H_
#include <AMReX_Config.H>

#include <AMReX_TypeTraits.H>
#include <AMReX_MultiFab.H>

namespace amrex
{

template <class PC, class MF, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
void
ParticleToMesh (PC const& pc, MF& mf, int lev, F&& f, bool zero_out_input=true)
{
    BL_PROFILE("amrex::ParticleToMesh");

    if (zero_out_input) { mf.setVal(0.0); }

    MultiFab mf_tmp(pc.ParticleBoxArray(lev),
                    pc.ParticleDistributionMap(lev),
                    mf.nComp(), mf.nGrowVect());
    mf_tmp.setVal(0.0);

    using ParIter = typename PC::ParConstIterType;
    const auto& plevel = pc.GetParticles(lev);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        for(ParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            const auto& tile = pti.GetParticleTile();
            const auto np = tile.numParticles();
            const auto& aos = tile.GetArrayOfStructs();
            const auto pstruct = aos().dataPtr();

            FArrayBox& fab = mf_tmp[pti];
            auto fabarr = fab.array();

            AMREX_FOR_1D( np, i,
            {
                f(pstruct[i], fabarr);
            });
        }
    }
    else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            FArrayBox local_fab;
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();

                FArrayBox& fab = mf_tmp[pti];

                Box tile_box = pti.tilebox();
                tile_box.grow(mf_tmp.nGrowVect());
                local_fab.resize(tile_box,mf_tmp.nComp());
                local_fab.setVal<RunOn::Host>(0.0);
                auto fabarr = local_fab.array();

                AMREX_FOR_1D( np, i,
                {
                    f(pstruct[i], fabarr);
                });

                fab.atomicAdd<RunOn::Host>(local_fab, tile_box, tile_box, 0, 0, mf_tmp.nComp());
            }
        }
    }

    mf.ParallelAdd(mf_tmp, 0, 0, mf_tmp.nComp(), mf_tmp.nGrowVect(), IntVect(0), pc.Geom(lev).periodicity());
}

template <class PC, class MF, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
void
MeshToParticle (PC& pc, MF const& mf, int lev, F&& f)
{
    BL_PROFILE("amrex::MeshToParticle");

    MultiFab* mf_pointer = pc.OnSameGrids(lev, mf) ?
        const_cast<MultiFab*>(&mf) : new MultiFab(pc.ParticleBoxArray(lev),
                                                  pc.ParticleDistributionMap(lev),
                                                  mf.nComp(), mf.nGrowVect());

    if (mf_pointer != &mf) mf_pointer->copy(mf,0,0,mf.nComp(),mf.nGrowVect(),mf.nGrowVect());

    auto& plevel = pc.GetParticles(lev);
    using ParIter = typename PC::ParIterType;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for(ParIter pti(pc, lev); pti.isValid(); ++pti)
    {
        auto& tile = pti.GetParticleTile();
        const auto np = tile.numParticles();
        auto& aos = tile.GetArrayOfStructs();
        auto pstruct = aos().dataPtr();

        const FArrayBox& fab = (*mf_pointer)[pti];
        auto fabarr = fab.array();

        AMREX_FOR_1D( np, i,
        {
            f(pstruct[i], fabarr);
        });
    }

    if (mf_pointer != &mf) delete mf_pointer;
}

}
#endif
