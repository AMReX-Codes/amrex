#ifndef AMREX_PARTICLEMESH_H_
#define AMREX_PARTICLEMESH_H_

namespace amrex
{

template <class PC, class MF, class F>
void
ParticleToMesh(PC const& pc, MF& mf, int lev, F f)
{
    BL_PROFILE("amrex::ParticleToMesh");
    
    MultiFab* mf_pointer = pc.OnSameGrids(lev, mf) ?
        &mf : new MultiFab(pc.ParticleBoxArray(lev), 
                           pc.ParticleDistributionMap(lev),
                           mf.nComp(), mf.nGrow());
    mf_pointer->setVal(0.);

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

            FArrayBox& fab = (*mf_pointer)[pti];
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
#ifdef _OPENMP
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

                FArrayBox& fab = (*mf_pointer)[pti];

                Box tile_box = pti.tilebox();
                tile_box.grow(mf_pointer->nGrow());
                local_fab.resize(tile_box,mf_pointer->nComp());
                local_fab = 0.0;
                auto fabarr = local_fab.array();
                
                AMREX_FOR_1D( np, i,
                {
                    f(pstruct[i], fabarr);
                });
                
                fab.atomicAdd(local_fab, tile_box, tile_box, 0, 0, mf_pointer->nComp());
            }
        }

        mf_pointer->SumBoundary(pc.Geom(lev).periodicity());

        if (mf_pointer != &mf)
        {
            mf.copy(*mf_pointer,0,0,mf_pointer->nComp());
            delete mf_pointer;
        }
    }
}

template <class PC, class MF, class F>
void
MeshToParticle(PC& pc, MF const& mf, int lev, F f)
{
    BL_PROFILE("amrex::MeshToParticle");

    MultiFab* mf_pointer = pc.OnSameGrids(lev, mf) ?
        const_cast<MultiFab*>(&mf) : new MultiFab(pc.ParticleBoxArray(lev), 
                                                  pc.ParticleDistributionMap(lev),
                                                  mf.nComp(), mf.nGrow());

    if (mf_pointer != &mf) mf_pointer->copy(mf,0,0,mf.nComp(),0,mf.nGrow());
        
    auto& plevel = pc.GetParticles(lev);
    using ParIter = typename PC::ParIterType;
#ifdef _OPENMP
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
