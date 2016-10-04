
#ifndef _ParticleContainer_H_
#define _ParticleContainer_H_

#include "Particles.H"

    // Particle Attributes:
    //    0 -- w  (weight)
    //    1 -- px (x mom)
    //    2 -- py (y mom)
    //    3 -- pz (z mom)
    //    4 -- Ex
    //    5 -- Ey
    //    6 -- Ez
    //    7 -- Bx
    //    8 -- By
    //    9 -- Bz
    //   10 -- gaminv


class MyParticleContainer
    : public ParticleContainer<11,0>
{

public:

    MyParticleContainer (const Array<Geometry>            & geom,
                         const Array<DistributionMapping> & dmap,
                         const Array<BoxArray>            & ba,
                         const Array<int>                 & rr)

     : ParticleContainer<11> (geom,dmap,ba,rr)
    { }

    virtual void Init(MultiFab& dummy_mf); 

    virtual void ChargeDeposition(MultiFab& mf_to_be_filled, int lev, int order) const;

    virtual void CurrentDeposition(PArray<MultiFab>& mf_to_be_filled, int lev,
                                   Real dt) const;

    virtual void FieldGather(MultiFab& Ex, MultiFab& Ey, MultiFab& Ez,
                             MultiFab& Bx, MultiFab& By, MultiFab& Bz,
                             long order, long field_gathe_algo);

    virtual void ParticlePush(Real dt);

};
#endif /*_ParticleContainer_H_*/
