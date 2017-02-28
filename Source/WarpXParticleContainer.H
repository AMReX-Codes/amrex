#ifndef WARPX_WarpXParticleContainer_H_
#define WARPX_WarpXParticleContainer_H_

#include <memory>

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

class MultiParticleContainer;

class WarpXParticleContainer
    : public ParticleContainer<PIdx::nattribs,0,std::vector<Particle<PIdx::nattribs,0> > >
{
public:
    friend MultiParticleContainer;

    WarpXParticleContainer (AmrCore* amr_core, int ispecies);
    virtual ~WarpXParticleContainer() {}

    virtual void AllocData () = 0;

    virtual void InitData () = 0;

    virtual void Evolve (int lev,
			 const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
			 const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
			 MultiFab& jx, MultiFab& jy, MultiFab& jz, Real t, Real dt) = 0;

    virtual void PushX (int lev, Real dt) = 0;

    virtual void PostRestart () = 0;

    std::unique_ptr<MultiFab> GetChargeDensity(int lev, bool local = false);

    void AddNParticles (int n, const Real* x, const Real* y, const Real* z,
			const Real* vx, const Real* vy, const Real* vz,
			int nattr, const Real* attr, int uniqueparticles);

    void ReadHeader (std::istream& is);

    void WriteHeader (std::ostream& os) const;

    const ParGDBBase& GDB () const { return *m_gdb; }

    static void ReadParameters ();

protected:

    int species_id;

    Real charge;
    Real mass;

    static bool do_tiling;
    static IntVect tile_size;
};

#endif
