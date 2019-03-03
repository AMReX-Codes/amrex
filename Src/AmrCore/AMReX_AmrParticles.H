#ifndef _AmrParticles_H_
#define _AmrParticles_H_

#include <AMReX_Particles.H>
#include <AMReX_TracerParticles.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>

namespace amrex {

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AssignDensity(int rho_index,
                                                                                 Vector<std::unique_ptr<MultiFab> >& mf_to_be_filled, 
                                                                                 int lev_min, int ncomp, int finest_level, int ngrow) const
{
    
    BL_PROFILE("ParticleContainer::AssignDensity()");
    
    if (rho_index != 0) amrex::Abort("AssignDensity only works if rho_index = 0");
    
    BL_ASSERT(NStructReal >= 1);
    BL_ASSERT(NStructReal >= ncomp);
    BL_ASSERT(ncomp == 1 || ncomp == AMREX_SPACEDIM+1);
    
    if (finest_level == -1) {
        finest_level = this->finestLevel();
    }
    while (!this->m_gdb->LevelDefined(finest_level)) {
        finest_level--;
    }

    ngrow = std::max(ngrow, 2);
    
    // Create the space for mf_to_be_filled, regardless of whether we'll need a temporary mf
    mf_to_be_filled.resize(finest_level+1);
    for (int lev = lev_min; lev <= finest_level; lev++)
    {        
        auto ng = lev == lev_min ? IntVect(AMREX_D_DECL(ngrow,ngrow,ngrow)) : this->m_gdb->refRatio(lev-1);
        mf_to_be_filled[lev].reset(new MultiFab(this->m_gdb->boxArray(lev),
                                                this->m_gdb->DistributionMap(lev), ncomp, ng));
        mf_to_be_filled[lev]->setVal(0.0);
    }
    
    // Test whether the grid structure of the boxArray is the same as the ParticleBoxArray at all levels 
    bool all_grids_the_same = true; 
    for (int lev = lev_min; lev <= finest_level; lev++)
    {
        if (!this->OnSameGrids(lev, *mf_to_be_filled[lev]))
        {
            all_grids_the_same = false;
            break;
        }
    }
    
    Vector<std::unique_ptr<MultiFab> > mf_part;
    if (!all_grids_the_same)
    { 
        // Create the space for the temporary, mf_part
        mf_part.resize(finest_level+1);
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            auto ng = lev == lev_min ? IntVect(AMREX_D_DECL(ngrow,ngrow,ngrow)) : this->m_gdb->refRatio(lev-1);
            mf_part[lev].reset(new MultiFab(this->ParticleBoxArray(lev), 
                                            this->ParticleDistributionMap(lev), ncomp, ng));
            mf_part[lev]->setVal(0.0);
        }
    }
    
    auto & mf = (all_grids_the_same) ? mf_to_be_filled : mf_part;
    
    if (finest_level == 0)
    {
        //
        // Just use the far simpler single-level version.
        //
        this->AssignCellDensitySingleLevel(rho_index, *mf[0], 0, ncomp);
        //
        // I believe that we don't need any information in ghost cells so we don't copy those.
        //
        if ( ! all_grids_the_same) {
            mf_to_be_filled[0]->copy(*mf[0],0,0,ncomp);
        }
        return;
    }

    // configure this to do a no-op at the physical boundaries.
    PhysBCFunct<BndryFuncArray> cphysbc, fphysbc;
    int lo_bc[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int hi_bc[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
    PCInterp mapper;
    
    Vector<std::unique_ptr<MultiFab> > tmp(finest_level+1);
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        const BoxArray& ba = mf[lev]->boxArray();
        const DistributionMapping& dm = mf[lev]->DistributionMap();
        tmp[lev].reset(new MultiFab(ba, dm, 1, 0));
        tmp[lev]->setVal(0.0);
    }
    
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        
        this->AssignCellDensitySingleLevel(rho_index, *mf[lev], lev, 1, 0);
        
        if (lev < finest_level) {
            PhysBCFunct<BndryFuncArray> cphysbc(this->m_gdb->Geom(lev), bcs,
                                                BndryFuncArray([](Real* data,
                                                                  AMREX_ARLIM_P(lo), AMREX_ARLIM_P(hi),
                                                                  const int* dom_lo, const int* dom_hi,
                                                                  const Real* dx, const Real* grd_lo,
                                                                  const Real* time, const int* bc){}));
            PhysBCFunct<BndryFuncArray> fphysbc(this->m_gdb->Geom(lev+1), bcs,
                                                BndryFuncArray([](Real* data,
                                                                  AMREX_ARLIM_P(lo), AMREX_ARLIM_P(hi),
                                                                  const int* dom_lo, const int* dom_hi,
                                                                  const Real* dx, const Real* grd_lo,
                                                                  const Real* time, const int* bc){}));

            amrex::InterpFromCoarseLevel(*tmp[lev+1], 0.0, *mf[lev],
                                         rho_index, rho_index, ncomp,
                                         this->m_gdb->Geom(lev),
                                         this->m_gdb->Geom(lev+1),
                                         cphysbc, 0, fphysbc, 0,
                                         this->m_gdb->refRatio(lev), &mapper,
                                         bcs, 0);
        }
        
        if (lev > lev_min) {
            // Note - this will double count the mass on the coarse level in 
            // regions covered by the fine level, but this will be corrected
            // below in the call to average_down.
            amrex::sum_fine_to_coarse(*mf[lev], *mf[lev-1], rho_index, 1,
                                      this->m_gdb->refRatio(lev-1),
                                      this->m_gdb->Geom(lev-1),
                                      this->m_gdb->Geom(lev));
        }
        
        mf[lev]->plus(*tmp[lev], rho_index, ncomp, 0);
    }
    
    for (int lev = finest_level - 1; lev >= lev_min; --lev) {
        amrex::average_down(*mf[lev+1], *mf[lev], rho_index, ncomp, this->m_gdb->refRatio(lev));
    }
    
    if (!all_grids_the_same) {
        for (int lev = lev_min; lev <= finest_level; lev++) {
            mf_to_be_filled[lev]->copy(*mf_part[lev],0,0,1);
        }
    }
    if (lev_min > 0) {
        int nlevels = finest_level - lev_min + 1;
        for (int i = 0; i < nlevels; i++)
            {
                mf_to_be_filled[i] = std::move(mf_to_be_filled[i+lev_min]);
            }
        mf_to_be_filled.resize(nlevels);
    }
}

template <int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0>
class AmrParticleContainer
        : public ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
{

public:

    typedef Particle<NStructReal, NStructInt> ParticleType;
    
    AmrParticleContainer (AmrCore* amr_core)
        : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>(amr_core->GetParGDB())
    {
    }

    AmrParticleContainer (const Vector<Geometry>            & geom, 
                          const Vector<DistributionMapping> & dmap,
                          const Vector<BoxArray>            & ba,
                          const Vector<int>                 & rr)
        : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>(geom, dmap, ba, rr)
    {
    }
    
    ~AmrParticleContainer () {}
    
};

class AmrTracerParticleContainer
    : public TracerParticleContainer
{
public:

    AmrTracerParticleContainer (AmrCore* amr_core)
	: TracerParticleContainer(amr_core->GetParGDB())
    {
    }

    ~AmrTracerParticleContainer () {}
};

}

#endif
