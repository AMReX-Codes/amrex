#ifndef _AmrParticles_H_
#define _AmrParticles_H_

#include <Particles.H>
#include <TracerParticles.H>
#include <AmrParGDB.H>

template <int NR, int NI=0>
class AmrParticleContainer
    : public ParticleContainer<NR,NI>
{

public:

    typedef Particle<NR,NI> ParticleType;
    typedef typename ParticleContainer<NR,NI>::PBox PBox;
    typedef typename ParticleContainer<NR,NI>::PMap PMap;
    
    AmrParticleContainer (Amr* amr)
	: ParticleContainer<NR,NI>(amr->GetParGDB())
    {
	this->m_particles.resize(this->m_gdb->finestLevel()+1);
    }

    ~AmrParticleContainer () {}
};

class AmrTracerParticleContainer
    : public TracerParticleContainer
{
public:

    AmrTracerParticleContainer (Amr* amr)
	: TracerParticleContainer(amr->GetParGDB())
    {
	this->m_particles.resize(this->m_gdb->finestLevel()+1);
    }

    ~AmrTracerParticleContainer () {}
};

#endif
