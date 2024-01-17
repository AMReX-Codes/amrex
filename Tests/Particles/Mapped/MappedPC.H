#ifndef MAPPED_PC_H_
#define MAPPED_PC_H_

#include <AMReX_Particles.H>

struct MappedRealIdx
{
    enum {
        vx = 0,
        vy, vz,
        ncomps
    };
};

struct MappedIntIdx
{
    enum {
        i = 0,
        j, k,
        ncomps
    };
};

struct MappedAssignor
{
    template <typename P>
    AMREX_GPU_HOST_DEVICE
    amrex::IntVect operator() (P const& p,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& /*plo*/,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& /*dxi*/,
                               const amrex::Box& /*domain*/) const noexcept
    {
        amrex::IntVect iv( AMREX_D_DECL(p.idata(0),
                                        p.idata(1),
                                        p.idata(2)));
        return iv;
    }
};

class MappedPC
    : public amrex::ParticleContainer<MappedRealIdx::ncomps, MappedIntIdx::ncomps, 0, 0,
                                      amrex::DefaultAllocator, MappedAssignor>
{

public:

    MappedPC (amrex::ParGDBBase* gdb)
        : amrex::ParticleContainer<MappedRealIdx::ncomps, MappedIntIdx::ncomps, 0, 0,
                                   amrex::DefaultAllocator, MappedAssignor>(gdb)
        {}

    MappedPC (const amrex::Geometry            & geom,
              const amrex::DistributionMapping & dmap,
              const amrex::BoxArray            & ba)
        : amrex::ParticleContainer<MappedRealIdx::ncomps, MappedIntIdx::ncomps, 0, 0,
                                   amrex::DefaultAllocator, MappedAssignor>(geom, dmap, ba)
        {}

    void InitParticles (amrex::MultiFab& a_xyz_loc);

    void AdvectWithUmac (amrex::MultiFab* umac, int level, amrex::Real dt,
                         const amrex::MultiFab& a_xyz_loc);
    void AdvectWithUCC  (amrex::MultiFab& u_cc, int level, amrex::Real dt,
                         const amrex::MultiFab& a_xyz_loc);
    void AdvectWithUND  (amrex::MultiFab& u_cc, int level, amrex::Real dt,
                         const amrex::MultiFab& a_xyz_loc);
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
