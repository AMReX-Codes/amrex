
#ifndef _ParticleContainer_H_
#define _ParticleContainer_H_

#include "Particles.H"

struct PIdx
{
    enum {
	// Particle Attributes:
	w = 0,  // weight
	px, py, pz, Ex, Ey, Ez, Bx, By, Bz, gaminv,
	nattribs
    };
};

class MyParticleContainer
    : public ParticleContainer<PIdx::nattribs,0>
{

public:

    MyParticleContainer (const Array<Geometry>            & geom,
                         const Array<DistributionMapping> & dmap,
                         const Array<BoxArray>            & ba,
                         const Array<int>                 & rr)

	: ParticleContainer<PIdx::nattribs,0> (geom,dmap,ba,rr)
    { }

    virtual ~MyParticleContainer() {}

    void Init(MultiFab& dummy_mf); 

    void ChargeDeposition(MultiFab& mf_to_be_filled, int lev, int order) const;

    void CurrentDeposition(PArray<MultiFab>& mf_to_be_filled, int lev,
			   Real dt) const;

    void FieldGather(MultiFab& Ex, MultiFab& Ey, MultiFab& Ez,
		     MultiFab& Bx, MultiFab& By, MultiFab& Bz,
		     long order, long field_gathe_algo);

    void ParticlePush(Real dt);
};
#endif /*_ParticleContainer_H_*/
