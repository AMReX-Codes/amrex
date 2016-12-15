#ifndef _AmrParticles_H_
#define _AmrParticles_H_

#include <Particles.H>
#include <TracerParticles.H>
#include <AmrParGDB.H>

template <int NR, int NI=0, class C = std::deque<Particle<NR,NI> > >
class AmrParticleContainer
    : public ParticleContainer<NR,NI,C>
{

public:

    typedef Particle<NR,NI> ParticleType;
    typedef typename ParticleContainer<NR,NI,C>::PBox PBox;
    typedef typename ParticleContainer<NR,NI,C>::PMap PMap;
    
    AmrParticleContainer (AmrCore* amr_core)
	: ParticleContainer<NR,NI,C>(amr_core->GetParGDB())
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

#endif
