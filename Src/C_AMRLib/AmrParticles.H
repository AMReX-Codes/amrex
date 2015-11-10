#ifndef _AmrParticles_H_
#define _AmrParticles_H_

#include <Particles.H>
#include <AmrParGDB.H>

template <int N>
class AmrParticleContainer
    : public ParticleContainer<N>
{

public:

    typedef Particle<N> ParticleType;
    typedef typename ParticleContainer<N>::PBox PBox;
    typedef typename ParticleContainer<N>::PMap PMap;
    
    AmrParticleContainer (Amr* amr)
	: ParticleContainer<N>(amr->GetParGDB())
    {
	this->m_particles.resize(this->m_gdb->finestLevel()+1);
    }

};

#endif
