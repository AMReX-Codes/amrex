#ifndef AMREX_AmrParticles_H_
#define AMREX_AmrParticles_H_
#include <AMReX_Config.H>

#include <AMReX_Particles.H>
#include <AMReX_TracerParticles.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>

namespace amrex {

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::AssignDensity (int rho_index,
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
        mf_to_be_filled[lev] = std::make_unique<MultiFab>(this->m_gdb->boxArray(lev),
                                                          this->m_gdb->DistributionMap(lev),
                                                          ncomp, ng);
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
            mf_part[lev] = std::make_unique<MultiFab>(this->ParticleBoxArray(lev),
                                                      this->ParticleDistributionMap(lev),
                                                      ncomp, ng);
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
        if ( ! all_grids_the_same) {
            mf_to_be_filled[0]->ParallelCopy(*mf[0],0,0,ncomp,0,0);
        }
        return;
    }

    // configure this to do a no-op at the physical boundaries.
    int lo_bc[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int hi_bc[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    Vector<BCRec> bcs(ncomp, BCRec(lo_bc, hi_bc));
    PCInterp mapper;

    Vector<std::unique_ptr<MultiFab> > tmp(finest_level+1);
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        const BoxArray& ba = mf[lev]->boxArray();
        const DistributionMapping& dm = mf[lev]->DistributionMap();
        tmp[lev] = std::make_unique<MultiFab>(ba, dm, ncomp, 0);
        tmp[lev]->setVal(0.0);
    }

    for (int lev = lev_min; lev <= finest_level; ++lev) {

        this->AssignCellDensitySingleLevel(rho_index, *mf[lev], lev, ncomp, 0);

        if (lev < finest_level) {
            PhysBCFunctNoOp cphysbc, fphysbc;
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
            amrex::sum_fine_to_coarse(*mf[lev], *mf[lev-1], rho_index, ncomp,
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
        //
        // We haven't interpolated the ghost cells so we can't copy them
        //
            mf_to_be_filled[lev]->ParallelCopy(*mf_part[lev],0,0,ncomp,0,0);
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

template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
void
ParticleToMesh (PC const& pc, const Vector<MultiFab*>& mf,
                int lev_min, int lev_max, F&& f,
                bool zero_out_input=true, bool vol_weight=true)
{
    BL_PROFILE("amrex::ParticleToMesh");

    if (lev_max == -1) { lev_max = pc.finestLevel(); }
    while (!pc.GetParGDB()->LevelDefined(lev_max)) { lev_max--; }

    if (lev_max == 0)
    {
        ParticleToMesh(pc, *mf[0], 0, std::forward<F>(f), zero_out_input);
        if (vol_weight) {
            const Real* dx = pc.Geom(0).CellSize();
            const Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
            mf[0]->mult(Real(1.0)/vol, 0, mf[0]->nComp(), mf[0]->nGrow());
        }
        return;
    }

    int ngrow = amrex::max(AMREX_D_DECL(mf[0]->nGrow(0),
                                        mf[0]->nGrow(1),
                                        mf[0]->nGrow(2)), 2);

    if (zero_out_input)
    {
        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            mf[lev]->setVal(Real(0.0));
        }
    }

    Vector<MultiFab> mf_part(lev_max+1);
    Vector<MultiFab> mf_tmp( lev_max+1);
    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
        mf_part[lev].define(pc.ParticleBoxArray(lev),
                            pc.ParticleDistributionMap(lev),
                            mf[lev]->nComp(), ngrow);
        mf_tmp[lev].define( pc.ParticleBoxArray(lev),
                            pc.ParticleDistributionMap(lev),
                            mf[lev]->nComp(), 0);
        mf_part[lev].setVal(0.0);
        mf_tmp[lev].setVal( 0.0);
    }

    // configure this to do a no-op at the physical boundaries.
    int lo_bc[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int hi_bc[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    Vector<BCRec> bcs(mf_part[lev_min].nComp(), BCRec(lo_bc, hi_bc));
    PCInterp mapper;

    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
        ParticleToMesh(pc, mf_part[lev], lev, std::forward<F>(f), zero_out_input);
        if (vol_weight) {
            const Real* dx = pc.Geom(lev).CellSize();
            const Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
            mf_part[lev].mult(Real(1.0)/vol, 0, mf_part[lev].nComp(), mf_part[lev].nGrow());
        }

        if (lev < lev_max) {
            PhysBCFunctNoOp cphysbc, fphysbc;
            amrex::InterpFromCoarseLevel(mf_tmp[lev+1], 0.0, mf_part[lev],
                                         0, 0, mf_part[lev].nComp(),
                                         pc.GetParGDB()->Geom(lev),
                                         pc.GetParGDB()->Geom(lev+1),
                                         cphysbc, 0, fphysbc, 0,
                                         pc.GetParGDB()->refRatio(lev), &mapper,
                                         bcs, 0);
        }

        if (lev > lev_min) {
            // Note - this will double count the mass on the coarse level in
            // regions covered by the fine level, but this will be corrected
            // below in the call to average_down.
            amrex::sum_fine_to_coarse(mf_part[lev], mf_part[lev-1], 0, mf_part[lev].nComp(),
                                      pc.GetParGDB()->refRatio(lev-1),
                                      pc.GetParGDB()->Geom(lev-1),
                                      pc.GetParGDB()->Geom(lev));
        }

        mf_part[lev].plus(mf_tmp[lev], 0, mf_part[lev].nComp(), 0);
    }

    for (int lev = lev_max-1; lev >= lev_min; --lev) {
        amrex::average_down(mf_part[lev+1], mf_part[lev], 0, mf_part[lev].nComp(),
                            pc.GetParGDB()->refRatio(lev));
    }

    for (int lev = lev_min; lev <= lev_max; lev++)
    {
        mf[lev]->ParallelCopy(mf_part[lev],0,0,mf_part[lev].nComp(),0,0);
    }
}

template <int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0,
          template<class> class Allocator=DefaultAllocator>
class AmrParticleContainer // NOLINT
    : public ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
{

public:

    typedef Particle<NStructReal, NStructInt> ParticleType;

    AmrParticleContainer () = default;

    AmrParticleContainer (AmrCore* amr_core)
        : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>(amr_core->GetParGDB())
    {
    }

    AmrParticleContainer (const Vector<Geometry>            & geom,
                          const Vector<DistributionMapping> & dmap,
                          const Vector<BoxArray>            & ba,
                          const Vector<int>                 & rr)
        : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>(geom, dmap, ba, rr)
    {
    }

    ~AmrParticleContainer () override = default;

    AmrParticleContainer ( const AmrParticleContainer &) = delete;
    AmrParticleContainer& operator= ( const AmrParticleContainer & ) = delete;

    AmrParticleContainer ( AmrParticleContainer && ) noexcept = default;
    AmrParticleContainer& operator= ( AmrParticleContainer && ) noexcept = default;
};

class AmrTracerParticleContainer
    : public TracerParticleContainer
{
public:

    AmrTracerParticleContainer (AmrCore* amr_core)
        : TracerParticleContainer(amr_core->GetParGDB())
    {
    }

    ~AmrTracerParticleContainer () override = default;

    AmrTracerParticleContainer ( const AmrTracerParticleContainer &) = delete;
    AmrTracerParticleContainer& operator= ( const AmrTracerParticleContainer & ) = delete;

    AmrTracerParticleContainer ( AmrTracerParticleContainer && ) noexcept = default;
    AmrTracerParticleContainer& operator= ( AmrTracerParticleContainer && ) noexcept = default;
};

}

#endif
