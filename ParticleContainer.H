
#ifndef _ParticleContainer_H_
#define _ParticleContainer_H_

#include <memory>
#include <map>

#include <Particles.H>
#include <AmrParGDB.H>

namespace{
  constexpr Real q_e = 1.6021766208e-19;
  constexpr Real m_e = 9.10938356e-31;
  constexpr Real m_p = 1.672621898e-27;
}

struct PIdx
{
    enum { // Particle Attributes stored in p.m_data
	w = 0,  // weight
	ux, uy, uz,
	nattribs
    };

    enum { // Particle data stored in partdata
	Ex=0, Ey, Ez, Bx, By, Bz,
	npartdata
    };
};

class MyParticleContainer
    : public ParticleContainer<PIdx::nattribs,0,std::vector<Particle<PIdx::nattribs,0> > >
{

public:

    MyParticleContainer (AmrCore* amr_core);

    virtual ~MyParticleContainer() {}

    void AllocData ();

    void InitData ();

    void Evolve (int lev,
		 const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
		 const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
		 MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt);

private:
    // The outermost Array is for AMR level.
    // map's key is grid id and value is an Array. The size of map is the size of local grids.
    // The Array (i.e., the value of the map) has npartdata unique_ptr to another Array.  
    // The size of the map value Array is npartdata.
    // The innermost Array holds a particular data (say Ex or gaminv) for particles on a grid.
    // The size of the innermost array is the number of particles.
    Array<std::map<int, Array<std::unique_ptr<Array<Real> > > > > m_partdata;

    Real charge;
    Real mass;
};
#endif /*_ParticleContainer_H_*/
