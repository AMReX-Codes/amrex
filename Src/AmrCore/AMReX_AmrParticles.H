#ifndef _AmrParticles_H_
#define _AmrParticles_H_

#include <AMReX_Particles.H>
#include <AMReX_TracerParticles.H>
#include <AMReX_AmrParGDB.H>

namespace amrex {

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
