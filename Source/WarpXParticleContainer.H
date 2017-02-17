#ifndef WARPX_WarpXParticleContainer_H_
#define WARPX_WarpXParticleContainer_H_

#include <memory>

#include <AMReX_Particles.H>
#include <AMReX_AmrParGDB.H>

struct PIdx
{
    enum { // Particle Attributes stored in amrex::ParticleContainer's struct of array
	w = 0,  // weight
	ux, uy, uz, Ex, Ey, Ez, Bx, By, Bz,
	nattribs
    };
};

class MultiParticleContainer;

class WarpXParticleContainer
    : public amrex::ParticleContainer<0,0,PIdx::nattribs>
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
			 amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz, amrex::Real t, amrex::Real dt) = 0;

    virtual void PostRestart () = 0;

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
};

#endif
