
#ifndef _ParticleContainer_H_
#define _ParticleContainer_H_

#include "Particles.H"

    // Particle Attributes:
    //    1 -- w  (weight)
    //    2 -- px (x mom)
    //    3 -- py (y mom)
    //    4 -- pz (z mom)
    //    5 -- Ex
    //    6 -- Ey
    //    7 -- Ez
    //    8 -- Bx
    //    9 -- By
    //   10 -- Bz


class MyParticleContainer
    : public ParticleContainer<10,0>
{

public:

    MyParticleContainer (const Array<Geometry>            & geom,
                         const Array<DistributionMapping> & dmap,
                         const Array<BoxArray>            & ba,
                         const Array<int>                 & rr)

     : ParticleContainer<10> (geom,dmap,ba,rr)
    { }

    virtual void ChargeDeposition(MultiFab& mf_to_be_filled, int lev, int order) const;

    virtual void CurrentDeposition(PArray<MultiFab>& mf_to_be_filled, int lev,
                                   Real dt) const;

    virtual void FieldGather(MultiFab& Ex, MultiFab& Ey, MultiFab& Ez,
                             MultiFab& Bx, MultiFab& By, MultiFab& Bz,
                             long order);

};
#endif /*_ParticleContainer_H_*/
