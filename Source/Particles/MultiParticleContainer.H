
#ifndef WARPX_ParticleContainer_H_
#define WARPX_ParticleContainer_H_

#include <memory>
#include <map>
#include <string>

#include <AMReX_Particles.H>

#include <WarpXParticleContainer.H>
#include <PhysicalParticleContainer.H>
#include <RigidInjectedParticleContainer.H>
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

    std::array<amrex::Real, 3> meanParticleVelocity(int ispecies) {
        return allcontainers[ispecies]->meanParticleVelocity();
    }

    void AllocData ();

    void InitData ();

#ifdef WARPX_DO_ELECTROSTATIC
    ///
    /// Performs the field gather operation using the input field E, for all the species
    /// in the MultiParticleContainer. This is the electrostatic version of the field gather.
    ///
    void FieldGatherES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                        const amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks);

    ///
    /// This evolves all the particles by one PIC time step, including charge deposition, the 
    /// field solve, and pushing the particles, for all the species in the MultiParticleContainer.
    /// This is the electrostatic version.
    ///
    void EvolveES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E, 
                         amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho, 
                   amrex::Real t, amrex::Real dt);

    ///
    /// This pushes the particle positions by one half time step for all the species in the 
    /// MultiParticleContainer. It is used to desynchronize the particles after initializaton
    /// or when restarting from a checkpoint. This is the electrostatic version.
    /// 
    void PushXES (amrex::Real dt);

    ///
    /// This deposits the particle charge onto rho, accumulating the value for all the species
    /// in the MultiParticleContainer. rho is assumed to contain node-centered multifabs. 
    /// This version is hard-coded for CIC deposition.
    ///
    void DepositCharge(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho, bool local = false);

    ///
    /// This returns the total particle charge for all the species in this MultiParticleContainer.
    /// This is needed to subtract the offset for periodic boundary conditions.
    ///
    amrex::Real sumParticleCharge(bool local = false);
#endif // WARPX_DO_ELECTROSTATIC
    
    ///
    /// Performs the field gather operation using the input fields E and B, for all the species
    /// in the MultiParticleContainer. This is the electromagnetic version of the field gather.
    ///
    void FieldGather (int lev,
                      const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
                      const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz); 

    ///
    /// This evolves all the particles by one PIC time step, including current deposition, the
    /// field solve, and pushing the particles, for all the species in the MultiParticleContainer.
    /// This is the electromagnetic version.
    ///
    void Evolve (int lev,
		 const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
		 const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
		 amrex::MultiFab& jx,  amrex::MultiFab& jy, amrex::MultiFab& jz, 
		 amrex::MultiFab* cjx,  amrex::MultiFab* cjy, amrex::MultiFab* cjz, 
                 amrex::MultiFab* rho, amrex::MultiFab* crho,
		 const amrex::MultiFab* cEx, const amrex::MultiFab* cEy, const amrex::MultiFab* cEz,
		 const amrex::MultiFab* cBx, const amrex::MultiFab* cBy, const amrex::MultiFab* cBz,
                 amrex::Real t, amrex::Real dt);

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
    /// This deposits the particle charge onto a node-centered MultiFab and returns a unique ptr
    /// to it. The charge density is accumulated over all the particles in the MultiParticleContainer
    /// This version uses PICSAR routines to perform the deposition.
    ///    
    std::unique_ptr<amrex::MultiFab> GetChargeDensity(int lev, bool local = false);

    void Checkpoint (const std::string& dir) const;

    void WritePlotFile( const std::string& dir,
                        const amrex::Vector<int>& real_flags, 
                        const amrex::Vector<std::string>& real_names) const;
    
    void Restart (const std::string& dir);

    void PostRestart ();

    void ReadHeader (std::istream& is);

    void WriteHeader (std::ostream& os) const;

    void SortParticlesByCell ();

    void Redistribute ();

    void RedistributeLocal (const int num_ghost);

    amrex::Vector<long> NumberOfParticlesInGrid(int lev) const;

    void Increment (amrex::MultiFab& mf, int lev);

    void SetParticleBoxArray (int lev, amrex::BoxArray& new_ba);
    void SetParticleDistributionMap (int lev, amrex::DistributionMapping& new_dm);

    int nSpecies() const {return nspecies;}

    int nSpeciesDepositOnMainGrid () const {
        int r = 0;
        for (int i : deposit_on_main_grid) {
            if (i) ++r;
        }
        return r;
    }

    void GetLabFrameData(const std::string& snapshot_name,
                         const int i_lab, const int direction,
                         const amrex::Real z_old, const amrex::Real z_new,
                         const amrex::Real t_boost, const amrex::Real t_lab, const amrex::Real dt,
                         amrex::Vector<WarpXParticleContainer::DiagnosticParticleData>& parts) const;
    
    // 
    // Parameters for the Cherenkov corrector in the FDTD solver.
    // Both stencils are calculated ar runtime.
    // 
    // Number of coefficients for the stencil of the NCI corrector.
    // The stencil is applied in the z direction only.
    static constexpr int nstencilz_fdtd_nci_corr=5;

    amrex::Vector<amrex::Array<amrex::Real, nstencilz_fdtd_nci_corr> > fdtd_nci_stencilz_ex;
    amrex::Vector<amrex::Array<amrex::Real, nstencilz_fdtd_nci_corr> > fdtd_nci_stencilz_by;

    std::vector<std::string> GetSpeciesNames() const { return species_names; }

    PhysicalParticleContainer& GetPCtmp () { return *pc_tmp; }
    
protected:

    // Particle container types
    enum struct PCTypes {Physical, RigidInjected};

    std::vector<std::string> species_names;

    std::vector<int> deposit_on_main_grid;

    std::vector<PCTypes> species_types;

private:

    // physical particles (+ laser)
    amrex::Vector<std::unique_ptr<WarpXParticleContainer> > allcontainers;
    // Temporary particle container, used e.g. for particle splitting.
    std::unique_ptr<PhysicalParticleContainer> pc_tmp;

    void ReadParameters ();

    // runtime parameters
    int nspecies = 1;   // physical particles only. If WarpX::use_laser, nspecies+1 == allcontainers.size().
};
#endif /*WARPX_ParticleContainer_H_*/
