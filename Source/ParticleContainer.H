
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

    ///
    /// Performs the field gather operation using the input field E, for all the species
    /// in the MultiParticleContainer. This is the electrostatic version of the field gather.
    ///
    void FieldGatherES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                        const amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks);

    ///
    /// Performs the field gather operation using the input fields E and B, for all the species
    /// in the MultiParticleContainer. This is the electromagnetic version of the field gather.
    ///
    void FieldGather (int lev,
                      const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
                      const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz); 

    ///
    /// This evolves all the particles by one PIC time step, including charge deposition, the 
    /// field solve, and pushing the particles, for all the species in the MultiParticleContainer.
    /// This is the electrostatic version.
    ///
    void EvolveES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E, 
                         amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho, 
                   amrex::Real t, amrex::Real dt);
   
    ///
    /// This evolves all the particles by one PIC time step, including current deposition, the
    /// field solve, and pushing the particles, for all the species in the MultiParticleContainer.
    /// This is the electromagnetic version.
    ///
    void Evolve (int lev,
		 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
		 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
		 amrex::MultiFab& jx,  amrex::MultiFab& jy, amrex::MultiFab& jz, 
                 amrex::MultiFab* rho, amrex::Real t, amrex::Real dt);

    ///
    /// This pushes the particle positions by one half time step for all the species in the 
    /// MultiParticleContainer. It is used to desynchronize the particles after initializaton
    /// or when restarting from a checkpoint. This is the electrostatic version.
    /// 
    void PushXES (amrex::Real dt);

    ///
    /// This pushes the particle positions by one half time step for all the species in the
    /// MultiParticleContainer. It is used to desynchronize the particles after initializaton
    /// or when restarting from a checkpoint. This is the electromagnetic version.
    /// 
    void PushX (amrex::Real dt);

    ///
    /// This pushes the particle momenta by dt for all the species in the
    /// MultiParticleContainer. It is used to desynchronize the particles after initializaton
    /// or when restarting from a checkpoint.  It is also used to synchronize particles at the
    /// the end of the run.  This is the electromagnetic version.
    /// 
    void PushP (int lev, amrex::Real dt,
                const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
                const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz);

    ///
    /// This deposits the particle charge onto rho, accumulating the value for all the species
    /// in the MultiParticleContainer. rho is assumed to contain node-centered multifabs. 
    /// This version is hard-coded for CIC deposition.
    ///
    void DepositCharge(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho, bool local = false);

    ///
    /// This deposits the particle charge onto a node-centered MultiFab and returns a unique ptr
    /// to it. The charge density is accumulated over all the particles in the MultiParticleContainer
    /// This version uses PICSAR routines to perform the deposition.
    ///    
    std::unique_ptr<amrex::MultiFab> GetChargeDensity(int lev, bool local = false);

    ///
    /// This returns the total particle charge for all the species in this MultiParticleContainer.
    /// This is needed to subtract the offset for periodic boundary conditions.
    ///
    amrex::Real sumParticleCharge(bool local = false);

    void Checkpoint (const std::string& dir,
		     const std::string& name,
		     bool is_checkpoint,
                     const amrex::Vector<std::string>& varnames = amrex::Vector<std::string>()) const;

    void Restart (const std::string& dir, const std::string& name);

    void PostRestart ();

    void ReadHeader (std::istream& is);

    void WriteHeader (std::ostream& os) const;

    void Redistribute ();

    void RedistributeLocal ();

    amrex::Vector<long> NumberOfParticlesInGrid(int lev) const;

    void Increment (amrex::MultiFab& mf, int lev);

    void SetParticleBoxArray (int lev, amrex::BoxArray& new_ba);
    void SetParticleDistributionMap (int lev, amrex::DistributionMapping& new_dm);

    int nSpecies() {return nspecies;}

protected:

    std::vector<std::string> species_names;

private:

    // physical particles (+ laser)
    amrex::Vector<std::unique_ptr<WarpXParticleContainer> > allcontainers;

    void ReadParameters ();

    // runtime parameters
    int nspecies = 1;   // physical particles only. If WarpX::use_laser, nspecies+1 == allcontainers.size().

};
#endif /*WARPX_ParticleContainer_H_*/
