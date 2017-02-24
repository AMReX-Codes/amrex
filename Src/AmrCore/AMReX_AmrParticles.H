#ifndef _AmrParticles_H_
#define _AmrParticles_H_

#include <AMReX_Particles.H>
#include <AMReX_TracerParticles.H>
#include <AMReX_AmrParGDB.H>

namespace amrex {

  template <int NR, int NI=0, int NA=0>
class AmrParticleContainer
    : public ParticleContainer<NR, NI, NA>
{

public:

    typedef Particle<NR,NI> ParticleType;
    
    AmrParticleContainer (AmrCore* amr_core)
      : ParticleContainer<NR, NI, NA>(amr_core->GetParGDB())
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
