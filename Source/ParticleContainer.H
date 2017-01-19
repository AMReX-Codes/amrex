
#ifndef _ParticleContainer_H_
#define _ParticleContainer_H_

#include <memory>
#include <map>

#include <Particles.H>
#include <AmrParGDB.H>

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

class MultiSpeciesContainer;

class SingleSpeciesContainer
    : public ParticleContainer<PIdx::nattribs,0,std::vector<Particle<PIdx::nattribs,0> > >
{
    friend MultiSpeciesContainer;

public:

    SingleSpeciesContainer (AmrCore* amr_core, int ispecies);
    virtual ~SingleSpeciesContainer() {}

    void AllocData ();

    void InitData ();

    void Evolve (int lev,
		 const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
		 const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
		 MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt);

    std::unique_ptr<MultiFab> GetChargeDensity(int lev, bool local = false);
 
    void AddNParticles (int n, const Real* x, const Real* y, const Real* z,
			const Real* vx, const Real* vy, const Real* vz,
			int nattr, const Real* attr, int uniqueparticles);
   
private:

    // The outermost Array is for AMR level.
    // map's key is grid id and value is an Array. The size of map is the size of local grids.
    // The Array (i.e., the value of the map) has npartdata unique_ptr to another Array.  
    // The size of the map value Array is npartdata.
    // The innermost Array holds a particular data (say Ex or gaminv) for particles on a grid.
    // The size of the innermost array is the number of particles.
    Array<std::map<int, Array<std::unique_ptr<Array<Real> > > > > m_partdata;

    int species_id;

    Real charge;
    Real mass;

    static bool do_tiling;
    static IntVect tile_size;
};

class MultiSpeciesContainer
{

public:

    MultiSpeciesContainer (AmrCore* amr_core);

    virtual ~MultiSpeciesContainer() {}

    SingleSpeciesContainer& GetSingleSpeciesContainer (int ispecies) {
	return *species[ispecies];
    }

    void AllocData ();

    void InitData ();

    void Evolve (int lev,
		 const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
		 const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
		 MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt);

    std::unique_ptr<MultiFab> GetChargeDensity(int lev, bool local = false);

    void Checkpoint (const std::string& dir, const std::string& name);

    void Restart (const std::string& dir, const std::string& name);

    void Redistribute (bool where_called=false);

private:

    Array<std::unique_ptr<SingleSpeciesContainer> > species;

    static void ReadStaticParameters ();

    // runtime parameters
    static bool do_tiling;
    static IntVect tile_size;
    static int nspecies;
};
#endif /*_ParticleContainer_H_*/
