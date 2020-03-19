/* Copyright 2019-2020 Andrew Myers, Ann Almgren, Axel Huebl
 * David Grote, Jean-Luc Vay, Junmin Gu
 * Luca Fedeli, Mathieu Lobet, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 * Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_ParticleContainer_H_
#define WARPX_ParticleContainer_H_

#include "WarpXParticleContainer.H"
#include "PhysicalParticleContainer.H"
#include "PhotonParticleContainer.H"
#include "Laser/LaserParticleContainer.H"
#include "Parser/WarpXParserWrapper.H"
#include "Particles/Collision/CollisionType.H"
#include "Particles/ParticleCreation/SmartUtils.H"
#include "Particles/ParticleCreation/SmartCopy.H"
#include "Particles/ParticleCreation/SmartCreate.H"
#include "Particles/ParticleCreation/FilterCopyTransform.H"
#include "Particles/RigidInjectedParticleContainer.H"

#ifdef WARPX_QED
#   include "Particles/ElementaryProcess/QEDInternals/QedChiFunctions.H"
#   include "Particles/ElementaryProcess/QEDInternals/BreitWheelerEngineWrapper.H"
#   include "Particles/ElementaryProcess/QEDInternals/QuantumSyncEngineWrapper.H"
#endif

#include <AMReX_Particles.H>

#include <memory>
#include <map>
#include <string>
#include <algorithm>

/**
 * The class MultiParticleContainer holds multiple instances of the polymorphic
 * class WarpXParticleContainer, stored in its member variable "allcontainers".
 * The class WarpX typically has a single (pointer to an) instance of
 * MultiParticleContainer.
 *
 * MultiParticleContainer typically has two types of functions:
 * - Functions that loop over all instances of WarpXParticleContainer in
 *   allcontainers and calls the corresponding function (for instance,
 *   MultiParticleContainer::Evolve loops over all particles containers and
 *   calls the corresponding WarpXParticleContainer::Evolve function).
 * - Functions that specifically handle multiple species (for instance
 *   ReadParameters or mapSpeciesProduct).
 */
class MultiParticleContainer
{

public:

    MultiParticleContainer (amrex::AmrCore* amr_core);

    ~MultiParticleContainer() {}

    WarpXParticleContainer& GetParticleContainer (int ispecies) {
        return *allcontainers[ispecies];
    }

#ifdef WARPX_USE_OPENPMD
    std::unique_ptr<WarpXParticleContainer>& GetUniqueContainer(int ispecies) {
      return  allcontainers[ispecies];
    }
#endif
    std::array<amrex::Real, 3> meanParticleVelocity(int ispecies) {
        return allcontainers[ispecies]->meanParticleVelocity();
    }

    void AllocData ();

    void InitData ();

    ///
    /// Performs the field gather operation using the input fields E and B, for all the species
    /// in the MultiParticleContainer. This is the electromagnetic version of the field gather.
    ///
    void FieldGather (int lev,
                      const amrex::MultiFab& Ex, const amrex::MultiFab& Ey,
                      const amrex::MultiFab& Ez, const amrex::MultiFab& Bx,
                      const amrex::MultiFab& By, const amrex::MultiFab& Bz);

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
                 amrex::Real t, amrex::Real dt, DtType a_dt_type=DtType::Full);

    ///
    /// This pushes the particle positions by one half time step for all the species in the
    /// MultiParticleContainer. It is used to desynchronize the particles after initializaton
    /// or when restarting from a checkpoint.
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
    ///
    std::unique_ptr<amrex::MultiFab> GetChargeDensity(int lev, bool local = false);

    void doFieldIonization ();

    void doCoulombCollisions ();

    void Checkpoint (const std::string& dir) const;

    void WritePlotFile (const std::string& dir) const;

    void Restart (const std::string& dir);

    void PostRestart ();

    void ReadHeader (std::istream& is);

    void WriteHeader (std::ostream& os) const;

    void SortParticlesByBin (amrex::IntVect bin_size);

    void Redistribute ();

    void RedistributeLocal (const int num_ghost);

    amrex::Vector<long> NumberOfParticlesInGrid(int lev) const;

    void Increment (amrex::MultiFab& mf, int lev);

    void SetParticleBoxArray (int lev, amrex::BoxArray& new_ba);
    void SetParticleDistributionMap (int lev, amrex::DistributionMapping& new_dm);

    int nSpecies() const {return nspecies;}

    int nSpeciesBackTransformedDiagnostics() const {return nspecies_back_transformed_diagnostics;}
    int mapSpeciesBackTransformedDiagnostics(int i) const {return map_species_back_transformed_diagnostics[i];}
    int doBackTransformedDiagnostics() const {return do_back_transformed_diagnostics;}

    int nSpeciesDepositOnMainGrid () const {
        bool const onMainGrid = true;
        auto const & v = m_deposit_on_main_grid;
        return std::count( v.begin(), v.end(), onMainGrid );
    }

    int nSpeciesGatherFromMainGrid() const {
        bool const fromMainGrid = true;
        auto const & v = m_gather_from_main_grid;
        return std::count( v.begin(), v.end(), fromMainGrid );
    }

    void GetLabFrameData(const std::string& snapshot_name,
                         const int i_lab, const int direction,
                         const amrex::Real z_old, const amrex::Real z_new,
                         const amrex::Real t_boost, const amrex::Real t_lab, const amrex::Real dt,
                         amrex::Vector<WarpXParticleContainer::DiagnosticParticleData>& parts) const;

    // Inject particles during the simulation (for particles entering the
    // simulation domain after some iterations, due to flowing plasma and/or
    // moving window).
    void ContinuousInjection(const amrex::RealBox& injection_box) const;
    // Update injection position for continuously-injected species.
    void UpdateContinuousInjectionPosition(amrex::Real dt) const;
    int doContinuousInjection() const;

    std::vector<std::string> GetSpeciesNames() const { return species_names; }

    PhysicalParticleContainer& GetPCtmp () { return *pc_tmp; }

    std::string m_B_ext_particle_s = "default";
    std::string m_E_ext_particle_s = "default";
    // External fields added to particle fields.
    amrex::Vector<amrex::Real> m_B_external_particle;
    amrex::Vector<amrex::Real> m_E_external_particle;
    // ParserWrapper for B_external on the particle
    std::unique_ptr<ParserWrapper<4> > m_Bx_particle_parser;
    std::unique_ptr<ParserWrapper<4> > m_By_particle_parser;
    std::unique_ptr<ParserWrapper<4> > m_Bz_particle_parser;
    // ParserWrapper for E_external on the particle
    std::unique_ptr<ParserWrapper<4> > m_Ex_particle_parser;
    std::unique_ptr<ParserWrapper<4> > m_Ey_particle_parser;
    std::unique_ptr<ParserWrapper<4> > m_Ez_particle_parser;

#ifdef WARPX_QED
    /**
    * \brief Performs QED events (Breit-Wheeler process and photon emission)
    */
    void doQedEvents();
#endif



protected:

#ifdef WARPX_QED
    /**
    * \brief Performs Breit-Wheeler process for the species for which it is enabled
    */
    void doQedBreitWheeler();

    /**
    * \brief Performs QED photon emission for the species for which it is enabled
    */
    void doQedQuantumSync();
#endif

    // Particle container types
    enum struct PCTypes {Physical, RigidInjected, Photon};

    std::vector<std::string> species_names;

    std::vector<std::string> lasers_names;

    std::vector<std::string> collision_names;

    amrex::Vector<std::unique_ptr<CollisionType> > allcollisions;

    //! instead of depositing (current, charge) on the finest patch level, deposit to the coarsest grid
    std::vector<bool> m_deposit_on_main_grid;

    //! instead of gathering fields from the finest patch level, gather from the coarsest
    std::vector<bool> m_gather_from_main_grid;

    std::vector<PCTypes> species_types;

    template<typename ...Args>
    amrex::MFItInfo getMFItInfo (const WarpXParticleContainer& pc_src,
                                 Args const&... pc_dsts) const noexcept
    {
        amrex::MFItInfo info;

        MFItInfoCheckTiling(pc_src, pc_dsts...);

        if (pc_src.do_tiling && amrex::Gpu::notInLaunchRegion()) {
            info.EnableTiling(pc_src.tile_size);
        }

#ifdef _OPENMP
        info.SetDynamic(true);
#endif

        return info;
    }


#ifdef WARPX_QED
    // The QED engines
    std::shared_ptr<BreitWheelerEngine> m_shr_p_bw_engine;
    std::shared_ptr<QuantumSynchrotronEngine> m_shr_p_qs_engine;
    //_______________________________

    /**
     * Initialize QED engines and provides smart pointers
     * to species who need QED processes
     */
    void InitQED ();

    //Variables to store how many species need a QED process
    int m_nspecies_quantum_sync = 0;
    int m_nspecies_breit_wheeler = 0;
    //________

    amrex::ParticleReal m_quantum_sync_photon_creation_energy_threshold =
      static_cast<amrex::ParticleReal>(
        2.0 * PhysConst::m_e * PhysConst::c * PhysConst::c );  /*!< Energy threshold for photon creation in Quantum Synchrotron process.*/

    /**
     * Returns the number of species having Quantum Synchrotron process enabled
     */
    int NSpeciesQuantumSync() const { return  m_nspecies_quantum_sync;}

    /**
     * Returns the number of species having Breit Wheeler process enabled
     */
    int NSpeciesBreitWheeler() const { return m_nspecies_breit_wheeler;}

    /**
     * Initializes the Quantum Synchrotron engine
     */
    void InitQuantumSync ();

    /**
     * Initializes the Quantum Synchrotron engine
     */
    void InitBreitWheeler ();

    /**
     * Called by InitQuantumSync if a new table has
     * to be generated.
     */
    void QuantumSyncGenerateTable();

    /**
     * Called by InitBreitWheeler if a new table has
     * to be generated.
     */
    void BreitWheelerGenerateTable();

#endif

private:

    // physical particles (+ laser)
    amrex::Vector<std::unique_ptr<WarpXParticleContainer> > allcontainers;
    // Temporary particle container, used e.g. for particle splitting.
    std::unique_ptr<PhysicalParticleContainer> pc_tmp;

    void ReadParameters ();

    void mapSpeciesProduct ();
    int getSpeciesID (std::string product_str);

    // Number of species dumped in BackTransformedDiagnostics
    int nspecies_back_transformed_diagnostics = 0;
    // map_species_back_transformed_diagnostics[i] is the species ID in
    // MultiParticleContainer for 0<i<nspecies_back_transformed_diagnostics
    std::vector<int> map_species_back_transformed_diagnostics;
    int do_back_transformed_diagnostics = 0;

    // runtime parameters
    int nlasers = 0;
    int nspecies = 1;   // physical particles only. nspecies+nlasers == allcontainers.size().
    int ncollisions = 0;

    void MFItInfoCheckTiling(const WarpXParticleContainer& pc_src) const noexcept
    {
        return;
    }

    template<typename First, typename ...Args>
    void MFItInfoCheckTiling(const WarpXParticleContainer& pc_src,
                             First const& pc_dst, Args const&... others) const noexcept
    {
        if (pc_src.do_tiling && amrex::Gpu::notInLaunchRegion()) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pc_dst.do_tiling,
                "For particle creation processes, either all or none of the "
                "particle species must use tiling.");
        }

        MFItInfoCheckTiling(pc_src, others...);
    }

    /**
    * Should be called right after mapSpeciesProduct in InitData.
    * It checks the physical correctness of product particle species
    * selected by the user for ionization process.
    */
    void CheckIonizationProductSpecies();

#ifdef WARPX_QED
    /**
    * Should be called right after mapSpeciesProduct in InitData.
    * It checks the physical correctness of product particle species
    * selected by the user for QED processes.
    */
    void CheckQEDProductSpecies();
#endif


};
#endif /*WARPX_ParticleContainer_H_*/
