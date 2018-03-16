#ifndef WARPX_WarpXParticleContainer_H_
#define WARPX_WarpXParticleContainer_H_

#include <memory>

#include <AMReX_Particles.H>
#include <AMReX_AmrCore.H>

struct PIdx
{
    enum { // Particle Attributes stored in amrex::ParticleContainer's struct of array
	w = 0,  // weight
	ux, uy, uz, Ex, Ey, Ez, Bx, By, Bz,
	nattribs
    };
};

class WarpXParIter
    : public amrex::ParIter<0,0,PIdx::nattribs>
{
public:
    using amrex::ParIter<0,0,PIdx::nattribs>::ParIter;

#if (BL_SPACEDIM == 2)
    void GetPosition (amrex::Vector<amrex::Real>& x,
                      amrex::Vector<amrex::Real>& y,
                      amrex::Vector<amrex::Real>& z) const;
    void SetPosition (const amrex::Vector<amrex::Real>& x,
                      const amrex::Vector<amrex::Real>& y,
                      const amrex::Vector<amrex::Real>& z);
#endif

    const std::array<amrex::Vector<amrex::Real>, PIdx::nattribs>& GetAttribs () const { 
        return GetStructOfArrays().GetRealData(); 
    }
    
    std::array<amrex::Vector<amrex::Real>, PIdx::nattribs>& GetAttribs () { 
        return GetStructOfArrays().GetRealData(); 
    }

    const amrex::Vector<amrex::Real>& GetAttribs (int comp) const { 
        return GetStructOfArrays().GetRealData(comp);
    }
    
    amrex::Vector<amrex::Real>& GetAttribs (int comp) { 
        return GetStructOfArrays().GetRealData(comp);
    }
};

class MultiParticleContainer;

class WarpXParticleContainer
    : public amrex::ParticleContainer<0,0,PIdx::nattribs>
{
public:
    friend MultiParticleContainer;

    WarpXParticleContainer (amrex::AmrCore* amr_core, int ispecies);
    virtual ~WarpXParticleContainer() {}

    virtual void InitData () = 0;

    virtual void FieldGatherES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                                const amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks) {}

    virtual void FieldGather (int lev,
                              const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
                              const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz) {}

    virtual void EvolveES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                                 amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho, 
                           amrex::Real t, amrex::Real dt) = 0;

    virtual void Evolve (int lev,
			 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
			 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
			 amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz,
                         amrex::MultiFab* rho, amrex::MultiFab* rho2,
                         amrex::Real t, amrex::Real dt) = 0;

    virtual void PostRestart () = 0;

    void AllocData ();

    ///
    /// This pushes the particle positions by one half time step.
    /// It is used to desynchronize the particles after initializaton
    /// or when restarting from a checkpoint.
    /// This is the electrostatic version of the particle push.
    /// 
    void PushXES (amrex::Real dt);

    ///
    /// This pushes the particle positions by one half time step.
    /// It is used to desynchronize the particles after initializaton
    /// or when restarting from a checkpoint.
    /// This is the electromagnetic version of the particle push.
    /// 
    void PushX (         amrex::Real dt);
    void PushX (int lev, amrex::Real dt);

    ///
    /// This pushes the particle momenta by dt.
    /// 
    virtual void PushP (int lev, amrex::Real dt,
                        const amrex::MultiFab& Ex,
                        const amrex::MultiFab& Ey,
                        const amrex::MultiFab& Ez,
                        const amrex::MultiFab& Bx,
                        const amrex::MultiFab& By,
                        const amrex::MultiFab& Bz) = 0;

    void DepositCharge(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                       bool local = false);
    std::unique_ptr<amrex::MultiFab> GetChargeDensity(int lev, bool local = false);

    /// 
    /// This returns the total charge for all the particles in this ParticleContainer.
    /// This is needed when solving Poisson's equation with periodic boundary conditions.
    ///
    amrex::Real sumParticleCharge(bool local = false);

    amrex::Real maxParticleVelocity(bool local = false);

    void AddNParticles (int lev,
                        int n, const amrex::Real* x, const amrex::Real* y, const amrex::Real* z,
			const amrex::Real* vx, const amrex::Real* vy, const amrex::Real* vz,
			int nattr, const amrex::Real* attr, int uniqueparticles);

    void AddOneParticle (int lev, int grid, int tile,
                         amrex::Real x, amrex::Real y, amrex::Real z,
                         const std::array<amrex::Real,PIdx::nattribs>& attribs);

    void AddOneParticle (ParticleTileType& particle_tile,
                         amrex::Real x, amrex::Real y, amrex::Real z,
                         const std::array<amrex::Real,PIdx::nattribs>& attribs);

    void ReadHeader (std::istream& is);

    void WriteHeader (std::ostream& os) const;

    static void ReadParameters ();

    static int NextID () { return ParticleType::NextID(); }

protected:

    int species_id;

    amrex::Real charge;
    amrex::Real mass;

    static int do_not_push;
};

#endif
