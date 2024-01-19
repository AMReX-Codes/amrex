#ifndef TERRAIN_FITTED_PC_H_
#define TERRAIN_FITTED_PC_H_

#include <AMReX_Particles.H>

struct TerrainRealIdx
{
    enum {
        vx = 0,
        vy, vz,
        ncomps
    };
};

struct TerrainIntIdx
{
    enum {
        k = 0,
        ncomps
    };
};

struct TerrainAssignor
{
    template <typename P>
    AMREX_GPU_HOST_DEVICE
    amrex::IntVect operator() (P const& p,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                               const amrex::Box& domain) const noexcept
    {
#if (AMREX_SPACEDIM == 2)
        amrex::IntVect iv( int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                           p.idata(0) );
        iv[0] += domain.smallEnd()[0];
#elif (AMREX_SPACEDIM == 3)
        amrex::IntVect iv( int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                           int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                           p.idata(0) );
        iv[0] += domain.smallEnd()[0];
        iv[1] += domain.smallEnd()[1];
#endif
        return iv;
    }
};

class TerrainPC
    : public amrex::ParticleContainer<TerrainRealIdx::ncomps, TerrainIntIdx::ncomps, 0, 0,
                                      amrex::DefaultAllocator, TerrainAssignor>
{

public:

    TerrainPC (amrex::ParGDBBase* gdb)
        : amrex::ParticleContainer<TerrainRealIdx::ncomps, TerrainIntIdx::ncomps, 0, 0,
                                   amrex::DefaultAllocator, TerrainAssignor>(gdb)
        {}

    TerrainPC (const amrex::Geometry            & geom,
               const amrex::DistributionMapping & dmap,
               const amrex::BoxArray            & ba)
        : amrex::ParticleContainer<TerrainRealIdx::ncomps, TerrainIntIdx::ncomps, 0, 0,
                                   amrex::DefaultAllocator, TerrainAssignor>(geom, dmap, ba)
        {}

    void InitParticles (amrex::MultiFab& a_z_loc);

    void AdvectWithUmac (amrex::MultiFab* umac, int level, amrex::Real dt, const amrex::MultiFab& a_z_loc);
    void AdvectWithUCC  (amrex::MultiFab& u_cc, int level, amrex::Real dt, const amrex::MultiFab& a_z_loc);
    void AdvectWithUND  (amrex::MultiFab& u_cc, int level, amrex::Real dt, const amrex::MultiFab& a_z_loc);

    void RedistributeLocal (bool remove_neg=true)
    {
        const int lev_min = 0;
        const int lev_max = finestLevel();
        const int nGrow = 0;
        const int local = 1;
        Redistribute(lev_min, lev_max, nGrow, local, remove_neg);
    }

    void RedistributeGlobal (bool remove_neg=true)
    {
        const int lev_min = 0;
        const int lev_max = finestLevel();
        const int nGrow = 0;
        const int local = 0;
        Redistribute(lev_min, lev_max, nGrow, local, remove_neg);
    }

};

#endif
