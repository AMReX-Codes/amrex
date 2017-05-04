
#ifndef WARPX_ParticleContainer_H_
#define WARPX_ParticleContainer_H_

#include <memory>
#include <map>

#include <AMReX_Particles.H>

#include <WarpXParticleContainer.H>
#include <PhysicalParticleContainer.H>
#include <LaserParticleContainer.H>

//
// MultiParticleContainer holds multiple (nspecies or npsecies+1 when
// using laser) WarpXParticleContainer.  WarpXParticleContainer is
// derived from amrex::ParticleContainer, and it is
// polymorphic. PhysicalParticleContainer and LaserParticleContainer are
// concrete class dervied from WarpXParticlContainer.
//

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

    void FieldGatherES (const amrex::Array<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E);

    void FieldGather (int lev,
                      const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
                      const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz); 

    ///
    /// This evolves all the particles by one PIC time step, including charge deposition, the 
    /// field solve, and pushing the particles. This is the electrostatic version.
    ///
    void EvolveES (const amrex::Array<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E, 
                         amrex::Array<std::unique_ptr<amrex::MultiFab> >& rho, 
                   amrex::Real t, amrex::Real dt);
   
    ///
    /// This evolves all the particles by one PIC time step, including current deposition, the
    /// field solve, and pushing the particles. This is the electromagnetic version.
    ///
    void Evolve (int lev,
		 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
		 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
		 amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz, 
                 amrex::Real t, amrex::Real dt);

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
    void PushX (int lev, amrex::Real dt);

    void DepositCharge(amrex::Array<std::unique_ptr<amrex::MultiFab> >& rho, bool local = false);
    std::unique_ptr<amrex::MultiFab> GetChargeDensity(int lev, bool local = false);

    amrex::Real sumParticleCharge(int lev, bool local = false);

    void Checkpoint (const std::string& dir,
		     const std::string& name,
		     bool is_checkpoint,
                     const amrex::Array<std::string>& varnames = amrex::Array<std::string>()) const;

    void Restart (const std::string& dir, const std::string& name);

    void PostRestart ();

    void ReadHeader (std::istream& is);

    void WriteHeader (std::ostream& os) const;

    void Redistribute ();

    amrex::Array<long> NumberOfParticlesInGrid(int lev) const;

    void Increment (amrex::MultiFab& mf, int lev);

    void SetParticleBoxArray (int lev, amrex::BoxArray& new_ba);
    void SetParticleDistributionMap (int lev, amrex::DistributionMapping& new_dm);

    int nSpecies() {return nspecies;}

protected:

    std::vector<std::string> species_names;

private:

    // physical particles (+ laser)
    amrex::Array<std::unique_ptr<WarpXParticleContainer> > allcontainers;

    void ReadParameters ();

    // runtime parameters
    int nspecies = 1;   // physical particles only. If WarpX::use_laser, nspecies+1 == allcontainers.size().

};
#endif /*WARPX_ParticleContainer_H_*/
