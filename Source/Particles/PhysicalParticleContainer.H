/* Copyright 2019-2020 Andrew Myers, Axel Huebl, David Grote
 * Ligia Diana Amorim, Luca Fedeli, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 * Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PhysicalParticleContainer_H_
#define WARPX_PhysicalParticleContainer_H_

#include "Initialization/PlasmaInjector.H"
#include "WarpXParticleContainer.H"
#include "Filter/NCIGodfreyFilter.H"
#include "Particles/ElementaryProcess/Ionization.H"

#include <AMReX_IArrayBox.H>

#ifdef WARPX_QED
    #include <QuantumSyncEngineWrapper.H>
    #include <BreitWheelerEngineWrapper.H>
    #include <QedChiFunctions.H>

    #include <Particles/ElementaryProcess/QEDPairGeneration.H>
    #include <Particles/ElementaryProcess/QEDPhotonEmission.H>
#endif

#include <map>

/**
 * PhysicalParticleContainer is the ParticleContainer class containing plasma
 * particles (if a simulation has 2 plasma species, say "electrons" and
 * "ions"), they will be two instances of PhysicalParticleContainer.
 *
 * PhysicalParticleContainer inherits from WarpXParticleContainer.
 */
class PhysicalParticleContainer
    : public WarpXParticleContainer
{
public:

    enum PhysicalParticleType{
        electron,
        positron,
        photon,
        other
    };

    PhysicalParticleContainer (amrex::AmrCore* amr_core,
                               int ispecies,
                               const std::string& name);

    PhysicalParticleContainer (amrex::AmrCore* amr_core);

    virtual ~PhysicalParticleContainer () {}

    virtual void InitData () override;

    void InitIonizationModule ();

    /**
     * \brief Apply external E and B fields on the particles. The E and B
     * fields could be defined as a constant or using a parser for reading
     * in a mathematical expression. The default value for the E- and B-fields
     * is (0.0,0.0,0.0).
     *
     * \param[in,out] Exp-Bzp pointer to fields on particles modified based
     *                  on external E and B
     * \param[in] xp,yp,zp arrays of particle positions required to compute
     *                  mathematical expression for the external fields
     *                  using parser.
     */
    void AssignExternalFieldOnParticles ( WarpXParIter& pti,
                                          RealVector& Exp, RealVector& Eyp,
                                          RealVector& Ezp, RealVector& Bxp,
                                          RealVector& Byp, RealVector& Bzp, int lev);

    virtual void FieldGather (int lev,
                              const amrex::MultiFab& Ex,
                              const amrex::MultiFab& Ey,
                              const amrex::MultiFab& Ez,
                              const amrex::MultiFab& Bx,
                              const amrex::MultiFab& By,
                              const amrex::MultiFab& Bz) final;

    void FieldGather (WarpXParIter& pti,
                      RealVector& Exp,
                      RealVector& Eyp,
                      RealVector& Ezp,
                      RealVector& Bxp,
                      RealVector& Byp,
                      RealVector& Bzp,
                      amrex::FArrayBox const * exfab,
                      amrex::FArrayBox const * eyfab,
                      amrex::FArrayBox const * ezfab,
                      amrex::FArrayBox const * bxfab,
                      amrex::FArrayBox const * byfab,
                      amrex::FArrayBox const * bzfab,
                      const int ngE, const int e_is_nodal,
                      const long offset,
                      const long np_to_gather,
                      int lev,
                      int depos_lev);

    /**
     * \brief Evolve is the central function PhysicalParticleContainer that
     * advances plasma particles for a time dt (typically one timestep).
     *
     * \param lev level on which particles are living
     * \param Ex MultiFab from which field Ex is gathered
     * \param Ey MultiFab from which field Ey is gathered
     * \param Ez MultiFab from which field Ez is gathered
     * \param Bx MultiFab from which field Bx is gathered
     * \param By MultiFab from which field By is gathered
     * \param Bz MultiFab from which field Bz is gathered
     * \param jx MultiFab to which the particles' current jx is deposited
     * \param jy MultiFab to which the particles' current jy is deposited
     * \param jz MultiFab to which the particles' current jz is deposited
     * \param cjx Same as jx (coarser, from lev-1), when using deposition buffers
     * \param cjy Same as jy (coarser, from lev-1), when using deposition buffers
     * \param cjz Same as jz (coarser, from lev-1), when using deposition buffers
     * \param rho MultiFab to which the particles' charge is deposited
     * \param crho Same as rho (coarser, from lev-1), when using deposition buffers
     * \param cEx Same as Ex (coarser, from lev-1), when using gather buffers
     * \param cEy Same as Ey (coarser, from lev-1), when using gather buffers
     * \param cEz Same as Ez (coarser, from lev-1), when using gather buffers
     * \param cBx Same as Bx (coarser, from lev-1), when using gather buffers
     * \param cBy Same as By (coarser, from lev-1), when using gather buffers
     * \param cBz Same as Bz (coarser, from lev-1), when using gather buffers
     * \param t current physical time
     * \param dt time step by which particles are advanced
     * \param a_dt_type type of time step (used for sub-cycling)
     *
     * Evolve iterates over particle iterator (each box) and performs filtering,
     * field gather, particle push and current deposition for all particles
     * in the box.
     */
    virtual void Evolve (int lev,
                         const amrex::MultiFab& Ex,
                         const amrex::MultiFab& Ey,
                         const amrex::MultiFab& Ez,
                         const amrex::MultiFab& Bx,
                         const amrex::MultiFab& By,
                         const amrex::MultiFab& Bz,
                         amrex::MultiFab& jx,
                         amrex::MultiFab& jy,
                         amrex::MultiFab& jz,
                         amrex::MultiFab* cjx,
                         amrex::MultiFab* cjy,
                         amrex::MultiFab* cjz,
                         amrex::MultiFab* rho,
                         amrex::MultiFab* crho,
                         const amrex::MultiFab* cEx,
                         const amrex::MultiFab* cEy,
                         const amrex::MultiFab* cEz,
                         const amrex::MultiFab* cBx,
                         const amrex::MultiFab* cBy,
                         const amrex::MultiFab* cBz,
                         amrex::Real t,
                         amrex::Real dt,
                         DtType a_dt_type=DtType::Full) override;

    virtual void PushPX (WarpXParIter& pti, amrex::Real dt, DtType a_dt_type=DtType::Full);

    virtual void PushP (int lev, amrex::Real dt,
                        const amrex::MultiFab& Ex,
                        const amrex::MultiFab& Ey,
                        const amrex::MultiFab& Ez,
                        const amrex::MultiFab& Bx,
                        const amrex::MultiFab& By,
                        const amrex::MultiFab& Bz) override;

    void PartitionParticlesInBuffers(
                        long& nfine_current,
                        long& nfine_gather,
                        long const np,
                        WarpXParIter& pti,
                        int const lev,
                        amrex::iMultiFab const* current_masks,
                        amrex::iMultiFab const* gather_masks,
                        RealVector& uxp,
                        RealVector& uyp,
                        RealVector& uzp,
                        RealVector& wp );

    void copy_attribs (WarpXParIter& pti);

    virtual void PostRestart () final {}

    void SplitParticles(int lev);

    IonizationFilterFunc getIonizationFunc ();

    // Inject particles in Box 'part_box'
    virtual void AddParticles (int lev);

    void AddPlasma(int lev, amrex::RealBox part_realbox = amrex::RealBox());

    void MapParticletoBoostedFrame(amrex::Real& x, amrex::Real& y, amrex::Real& z, std::array<amrex::Real, 3>& u);

    void AddGaussianBeam(amrex::Real x_m, amrex::Real y_m, amrex::Real z_m,
                         amrex::Real x_rms, amrex::Real y_rms, amrex::Real z_rms,
                         amrex::Real q_tot, long npart, int do_symmetrize);

    void CheckAndAddParticle(amrex::Real x, amrex::Real y, amrex::Real z,
                             std::array<amrex::Real, 3> u,
                             amrex::Real weight,
                             amrex::Gpu::HostVector<amrex::ParticleReal>& particle_x,
                             amrex::Gpu::HostVector<amrex::ParticleReal>& particle_y,
                             amrex::Gpu::HostVector<amrex::ParticleReal>& particle_z,
                             amrex::Gpu::HostVector<amrex::ParticleReal>& particle_ux,
                             amrex::Gpu::HostVector<amrex::ParticleReal>& particle_uy,
                             amrex::Gpu::HostVector<amrex::ParticleReal>& particle_uz,
                             amrex::Gpu::HostVector<amrex::ParticleReal>& particle_w);

    virtual void GetParticleSlice(const int direction, const amrex::Real z_old,
                                  const amrex::Real z_new, const amrex::Real t_boost,
                                  const amrex::Real t_lab, const amrex::Real dt,
                                  DiagnosticParticles& diagnostic_particles) final;

    virtual void ConvertUnits (ConvertDirection convert_dir) override;

/**
 * \brief Apply NCI Godfrey filter to all components of E and B before gather
 * \param lev MR level
 * \param box box onto which the filter is applied
 * \param exeli safeguard Elixir object (to avoid de-allocating too early
            --between ParIter iterations-- on GPU) for field Ex
 * \param eyeli safeguard Elixir object (to avoid de-allocating too early
            --between ParIter iterations-- on GPU) for field Ey
 * \param ezeli safeguard Elixir object (to avoid de-allocating too early
            --between ParIter iterations-- on GPU) for field Ez
 * \param bxeli safeguard Elixir object (to avoid de-allocating too early
            --between ParIter iterations-- on GPU) for field Bx
 * \param byeli safeguard Elixir object (to avoid de-allocating too early
            --between ParIter iterations-- on GPU) for field By
 * \param bzeli safeguard Elixir object (to avoid de-allocating too early
            --between ParIter iterations-- on GPU) for field Bz
 * \param filtered_Ex Array containing filtered value
 * \param filtered_Ey Array containing filtered value
 * \param filtered_Ez Array containing filtered value
 * \param filtered_Bx Array containing filtered value
 * \param filtered_By Array containing filtered value
 * \param filtered_Bz Array containing filtered value
 * \param Ex Field array before filtering (not modified)
 * \param Ey Field array before filtering (not modified)
 * \param Ez Field array before filtering (not modified)
 * \param Bx Field array before filtering (not modified)
 * \param By Field array before filtering (not modified)
 * \param Bz Field array before filtering (not modified)
 * \param exfab pointer to the Ex field (modified)
 * \param eyfab pointer to the Ey field (modified)
 * \param ezfab pointer to the Ez field (modified)
 * \param bxfab pointer to the Bx field (modified)
 * \param byfab pointer to the By field (modified)
 * \param bzfab pointer to the Bz field (modified)
 *
 * The NCI Godfrey filter is applied on Ex, the result is stored in filtered_Ex
 * and the pointer exfab is modified (before this function is called, it points to Ex
 * and after this function is called, it points to Ex_filtered)
 */
    void applyNCIFilter (
        int lev, const amrex::Box& box,
        amrex::Elixir& exeli, amrex::Elixir& eyeli, amrex::Elixir& ezeli,
        amrex::Elixir& bxeli, amrex::Elixir& byeli, amrex::Elixir& bzeli,
        amrex::FArrayBox& filtered_Ex, amrex::FArrayBox& filtered_Ey,
        amrex::FArrayBox& filtered_Ez, amrex::FArrayBox& filtered_Bx,
        amrex::FArrayBox& filtered_By, amrex::FArrayBox& filtered_Bz,
        const amrex::FArrayBox& Ex, const amrex::FArrayBox& Ey,
        const amrex::FArrayBox& Ez, const amrex::FArrayBox& Bx,
        const amrex::FArrayBox& By, const amrex::FArrayBox& Bz,
        amrex::FArrayBox const * & exfab, amrex::FArrayBox const * & eyfab,
        amrex::FArrayBox const * & ezfab, amrex::FArrayBox const * & bxfab,
        amrex::FArrayBox const * & byfab, amrex::FArrayBox const * & bzfab);

#ifdef WARPX_QED
    //Functions decleared in WarpXParticleContainer.H
    //containers for which QED processes could be relevant
    //are expected to override these functions

    /**
     * Tells if this PhysicalParticleContainer has Quantum
     * Synchrotron process enabled
     * @return true if process is enabled
     */
    bool has_quantum_sync() override;

    /**
     * Tells if this PhysicalParticleContainer has Breit
     * Wheeler process enabled
     * @return true if process is enabled
     */
    bool has_breit_wheeler() override;

    /**
     * Acquires a shared smart pointer to a BreitWheelerEngine
     * @param[in] ptr the pointer
     */
    void set_breit_wheeler_engine_ptr
        (std::shared_ptr<BreitWheelerEngine> ptr) override;

    /**
     * Acquires a shared smart pointer to a QuantumSynchrotronEngine
     * @param[in] ptr the pointer
     */
    void set_quantum_sync_engine_ptr
        (std::shared_ptr<QuantumSynchrotronEngine> ptr) override;
    //__________

    /**
     * This function evolves the optical depth of the particles if QED effects
     * are enabled.
     * @param[in,out] pti particle iterator (optical depth will be modified)
     * @param[in] dt temporal step
     */
    virtual void EvolveOpticalDepth(WarpXParIter& pti,
        amrex::Real dt);

    PhotonEmissionFilterFunc getPhotonEmissionFilterFunc ();

    PairGenerationFilterFunc getPairGenerationFilterFunc ();
#endif

protected:
    std::string species_name;
    std::unique_ptr<PlasmaInjector> plasma_injector;

    // When true, adjust the transverse particle positions accounting
    // for the difference between the Lorentz transformed time of the
    // particle and the time of the boosted frame.
    bool boost_adjust_transverse_positions = false;
    bool do_backward_propagation = false;

    // Inject particles during the whole simulation
    void ContinuousInjection (const amrex::RealBox& injection_box) override;

    //This function return true if the PhysicalParticleContainer contains electrons
    //or positrons, false otherwise

    //When true PhysicalParticleContainer tries to use a pusher including
    //radiation reaction
    bool do_classical_radiation_reaction = false;

#ifdef WARPX_QED
    // A flag to enable quantum_synchrotron process for leptons
    bool m_do_qed_quantum_sync = false;

    // A flag to enable breit_wheeler process [photons only!!]
    bool m_do_qed_breit_wheeler = false;

    // A smart pointer to an instance of a Quantum Synchrotron engine
    std::shared_ptr<QuantumSynchrotronEngine> m_shr_p_qs_engine;

    // A smart pointer to an instance of a Breit Wheeler engine [photons only!]
    std::shared_ptr<BreitWheelerEngine> m_shr_p_bw_engine;

#endif

};

#endif
