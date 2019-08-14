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

    const auto& plevel = pc.GetParticles(lev);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        for(MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());            
            if(plevel.find(index) == plevel.end()) continue;            
            const auto& tile = plevel.at(index);
            
            const auto np = tile.numParticles();
            if (np == 0) continue;
            
            const auto& aos = tile.GetArrayOfStructs();
            const auto pstruct = aos().dataPtr();
            
            FArrayBox& fab = (*mf_pointer)[mfi];
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
            for(MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
            {
                auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());            
                if(plevel.find(index) == plevel.end()) continue;            
                const auto& tile = plevel.at(index);
                
                const auto np = tile.numParticles();
                if (np == 0) continue;
                
                const auto& aos = tile.GetArrayOfStructs();
                const auto pstruct = aos().dataPtr();
                
                FArrayBox& fab = (*mf_pointer)[mfi];
                
                Box tile_box = mfi.tilebox();
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
    
}


}
#endif
