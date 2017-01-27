
#ifndef WARPX_ParticleContainer_H_
#define WARPX_ParticleContainer_H_

#include <memory>
#include <map>

#include <AMReX_Particles.H>
#include <AMReX_AmrParGDB.H>

//
// MultiParticleContainer holds multiple (nspecies or npsecies+1 when
// using laser) WarpXParticleContainer.  WarpXParticlContainer is
// derived from amrex::ParticleContainer, and it is
// polymorphic. PhysicalParticleContainer and LaserParticleContainer are
// concrete class dervied from WarpXParticlContainer.
//

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

class MultiParticleContainer;

class WarpXParticleContainer
    : public amrex::ParticleContainer<PIdx::nattribs,0,std::vector<amrex::Particle<PIdx::nattribs,0> > >
{
public:
    friend MultiParticleContainer;

    WarpXParticleContainer (amrex::AmrCore* amr_core, int ispecies);
    virtual ~WarpXParticleContainer() {}

    virtual void AllocData () = 0;

    virtual void InitData () = 0;

    virtual void Evolve (int lev,
			 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
			 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
			 amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz, amrex::Real dt) = 0;

    std::unique_ptr<amrex::MultiFab> GetChargeDensity(int lev, bool local = false);
 
    void AddNParticles (int n, const amrex::Real* x, const amrex::Real* y, const amrex::Real* z,
			const amrex::Real* vx, const amrex::Real* vy, const amrex::Real* vz,
			int nattr, const amrex::Real* attr, int uniqueparticles);

    void ReadHeader (std::istream& is);

    void WriteHeader (std::ostream& os) const;

    const amrex::ParGDBBase& GDB () const { return *m_gdb; }

    static void ReadParameters ();
   
protected:

    int species_id;

    amrex::Real charge;
    amrex::Real mass;

    static bool do_tiling;
    static amrex::IntVect tile_size;
};

class PhysicalParticleContainer
    : public WarpXParticleContainer
{
public:
    PhysicalParticleContainer (amrex::AmrCore* amr_core, int ispecies);
    virtual ~PhysicalParticleContainer () {}

    virtual void AllocData () override;
    
    virtual void InitData () override;

    virtual void Evolve (int lev,
			 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
			 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
			 amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz, amrex::Real dt) override;

private:
    // The outermost Array is for AMR level.
    // map's key is grid id and value is an Array. The size of map is the size of local grids.
    // The Array (i.e., the value of the map) has npartdata unique_ptr to another Array.  
    // The size of the map value Array is npartdata.
    // The innermost Array holds a particular data (say Ex or gaminv) for particles on a grid.
    // The size of the innermost array is the number of particles.
    amrex::Array<std::map<int, amrex::Array<std::unique_ptr<amrex::Array<amrex::Real> > > > > m_partdata;
};

class LaserParticleContainer
    : public WarpXParticleContainer
{
public:
    LaserParticleContainer (amrex::AmrCore* amr_core, int ispecies);
    virtual ~LaserParticleContainer () {}

    virtual void AllocData () override;
    
    virtual void InitData () override;

    virtual void Evolve (int lev,
			 const amrex::MultiFab&, const amrex::MultiFab&, const amrex::MultiFab&,
			 const amrex::MultiFab&, const amrex::MultiFab&, const amrex::MultiFab&,
			 amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz, amrex::Real dt) override;
};

class MultiParticleContainer
{

public:

    MultiParticleContainer (amrex::AmrCore* amr_core);

    ~MultiParticleContainer() {}

    WarpXParticleContainer& GetParticleContainer (int ispecies) {
	return *allcontainers[ispecies];
    }

    void AllocData ();

    void InitData ();

    void Evolve (int lev,
		 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
		 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
		 amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz, amrex::Real dt);

    std::unique_ptr<amrex::MultiFab> GetChargeDensity(int lev, bool local = false);

    void Checkpoint (const std::string& dir,
		     const std::string& name,
		     bool is_checkpoint) const;

    void Restart (const std::string& dir, const std::string& name);

    void ReadHeader (std::istream& is);

    void WriteHeader (std::ostream& os) const;

    void Redistribute (bool where_called=false);

private:

    // physical particles (+ laser)
    amrex::Array<std::unique_ptr<WarpXParticleContainer> > allcontainers;

    void ReadParameters ();

    // runtime parameters
    int nspecies = 1;   // physical particles only. If WarpX::use_laser, nspecies+1 == allcontainers.size().
};
#endif /*WARPX_ParticleContainer_H_*/
