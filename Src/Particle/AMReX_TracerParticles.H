#ifndef BL_TRACERPARTICLES_H_
#define BL_TRACERPARTICLES_H_

#include <AMReX_Particles.H>

namespace amrex {

class TracerParticleContainer
    : public ParticleContainer<BL_SPACEDIM>
{
public:

    TracerParticleContainer (ParGDBBase* gdb)
	: ParticleContainer<BL_SPACEDIM>(gdb)
	{}

    TracerParticleContainer (const Geometry            & geom, 
			     const DistributionMapping & dmap,
			     const BoxArray            & ba)
	: ParticleContainer<BL_SPACEDIM>(geom,dmap,ba)
	{}

    ~TracerParticleContainer () {}

    void AdvectWithUmac (MultiFab* umac, int level, Real dt);

    void AdvectWithUcc (const MultiFab& ucc, int level, Real dt);

    void Timestamp (const std::string& file, const MultiFab& mf, int lev, Real time,
		    const std::vector<int>& idx);
};

using TracerParIter = ParIter<BL_SPACEDIM>;

}

#endif

