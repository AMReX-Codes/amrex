
#ifndef _ParticleContainer_H_
#define _ParticleContainer_H_

#include <memory>
#include <map>

#include <Particles.H>

struct PIdx
{
    enum { // Particle Attributes stored in p.m_data
	w = 0,  // weight
	px, py, pz,
	nattribs
    };

    enum { // Particle data stored in partdata
	Ex=0, Ey, Ez, Bx, By, Bz, gaminv,
	npartdata
    };
};

class MyParticleContainer
    : public ParticleContainer<PIdx::nattribs,0>
{

public:

    MyParticleContainer (const Array<Geometry>            & geom,
                         const Array<DistributionMapping> & dmap,
                         const Array<BoxArray>            & ba,
                         const Array<int>                 & rr);

    virtual ~MyParticleContainer() {}

    void Init (MultiFab& dummy_mf); 

    void ChargeDeposition (MultiFab& mf_to_be_filled, int lev, int order) const;

    void CurrentDeposition (Array< std::unique_ptr<MultiFab> >& mf_to_be_filled,
			    int lev, Real dt) const;

    void GatherField (const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
		      const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz);

    void ParticlePush (Real dt);

private:
    // map's key is grid id and value is an Array. The size of map is the size of local grids.
    // Array has npartdata unique_ptr to another Array.  The size of the firs Array is npartdata.
    // The second Array holds a particular data (say Ex or gaminv) for particles on a grid.
    //      The size of the second array is the number of particles.
    std::map<int, Array<std::unique_ptr<Array<Real> > > > partdata;
};
#endif /*_ParticleContainer_H_*/
