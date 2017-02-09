#ifndef _AmrParticles_H_
#define _AmrParticles_H_

#include <AMReX_Particles.H>
#include <AMReX_TracerParticles.H>
#include <AMReX_AmrParGDB.H>

namespace amrex {

  template <int NR, int NI=0, int NCOMP=0>
class AmrParticleContainer
    : public ParticleContainer<NR, NI, NCOMP>
{

public:

  typedef Particle<NR,NI> ParticleType;
  typedef typename ParticleContainer<NR,NI,NCOMP>::PBox PBox;
  typedef typename ParticleContainer<NR,NI,NCOMP>::PMap PMap;
    
    AmrParticleContainer (AmrCore* amr_core)
      : ParticleContainer<NR, NI, NCOMP>(amr_core->GetParGDB())
    {
	this->m_particles.reserve(this->m_gdb->maxLevel()+1);	
	this->m_particles.resize(this->m_gdb->finestLevel()+1);
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
	this->m_particles.reserve(this->m_gdb->maxLevel()+1);
	this->m_particles.resize(this->m_gdb->finestLevel()+1);
    }

    ~AmrTracerParticleContainer () {}
};

}

#endif
