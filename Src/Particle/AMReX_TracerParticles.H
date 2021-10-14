#ifndef AMREX_TRACERPARTICLES_H_
#define AMREX_TRACERPARTICLES_H_
#include <AMReX_Config.H>

#include <AMReX_Particles.H>

namespace amrex {

class TracerParticleContainer
    : public ParticleContainer<AMREX_SPACEDIM>
{
public:

    TracerParticleContainer (ParGDBBase* gdb)
        : ParticleContainer<AMREX_SPACEDIM>(gdb)
        {}

    TracerParticleContainer (const Geometry            & geom,
                             const DistributionMapping & dmap,
                             const BoxArray            & ba)
        : ParticleContainer<AMREX_SPACEDIM>(geom,dmap,ba)
        {}

    ~TracerParticleContainer () {}

    TracerParticleContainer ( const TracerParticleContainer &) = delete;
    TracerParticleContainer& operator= ( const TracerParticleContainer & ) = delete;

    TracerParticleContainer ( TracerParticleContainer && ) = default;
    TracerParticleContainer& operator= ( TracerParticleContainer && ) = default;

    void AdvectWithUmac (MultiFab* umac, int level, Real dt);

    void AdvectWithUcc (const MultiFab& ucc, int level, Real dt);

    void Timestamp (const std::string& file, const MultiFab& mf, int lev, Real time,
                    const std::vector<int>& idx);
};

using TracerParIter = ParIter<AMREX_SPACEDIM>;

}

#endif

