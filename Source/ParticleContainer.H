
#ifndef _ParticleContainer_H_
#define _ParticleContainer_H_

#include <memory>
#include <map>

#include <AMReX_Particles.H>
#include <AMReX_AmrParGDB.H>

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
    : public amrex::ParticleContainer<PIdx::nattribs,0,std::vector<amrex::Particle<PIdx::nattribs,0> > >
{

public:

    MyParticleContainer (amrex::AmrCore* amr_core);

    virtual ~MyParticleContainer() {}

    void AllocData ();

    void InitData ();

    void Evolve (int lev,
		 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
		 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
		 amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz, amrex::Real dt);

    std::unique_ptr<amrex::MultiFab> GetChargeDensity(int lev, bool local = false);

private:

    // The outermost Array is for AMR level.
    // map's key is grid id and value is an Array. The size of map is the size of local grids.
    // The Array (i.e., the value of the map) has npartdata unique_ptr to another Array.  
    // The size of the map value Array is npartdata.
    // The innermost Array holds a particular data (say Ex or gaminv) for particles on a grid.
    // The size of the innermost array is the number of particles.
    amrex::Array<std::map<int, amrex::Array<std::unique_ptr<amrex::Array<amrex::Real> > > > > m_partdata;

    amrex::Real charge;
    amrex::Real mass;

    static void ReadStaticParameters ();

    // runtime parameters
    static int do_tiling;
    static amrex::IntVect tile_size;
};
#endif /*_ParticleContainer_H_*/
