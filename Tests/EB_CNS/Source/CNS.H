#ifndef CNS_H_
#define CNS_H_

#include <CNS_index_macros.H>
#include <CNS_parm.H>
#include <cns_prob_parm.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFluxRegister.H>

using namespace amrex;
class CNS
    :
    public amrex::AmrLevel
{
public:

    CNS ();
    CNS (amrex::Amr&            papa,
         int                    lev,
         const amrex::Geometry& level_geom,
         const amrex::BoxArray& bl,
         const amrex::DistributionMapping& dm,
         amrex::Real            time);
    virtual ~CNS ();

    CNS (const CNS& rhs) = delete;
    CNS& operator= (const CNS& rhs) = delete;

    // Restart from a checkpoint file.
    virtual void restart (amrex::Amr&     papa,
                          std::istream&   is,
                          bool            bReadSpecial = false) override;

    // Write checkpoint
    virtual void checkPoint(const std::string& dir,
                            std::ostream&      os,
                            amrex::VisMF::How  how = amrex::VisMF::NFiles,
                            bool               dump_old = true) override;

    virtual std::string thePlotFileType () const override {
        return {"HyperCLaw-V1.1"};
    }

    // Write a plotfile to specified directory.
    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                amrex::VisMF::How  how) override;

    // Initialize data on this level from another CNS (during regrid).
    virtual void init (amrex::AmrLevel& old) override;

    // Initialize data on this level after regridding if old level did not previously exist
    virtual void init () override;

    // Initialize grid data at problem start-up.
    virtual void initData () override;

    // Advance grids at this level in time.
    virtual amrex::Real advance (amrex::Real time,
                                 amrex::Real dt,
                                 int  iteration,
                                 int  ncycle) override;

    virtual void computeInitialDt (int                                 finest_level,
                                   int                                 sub_cycle,
                                   amrex::Vector<int>&                  n_cycle,
                                   const amrex::Vector<amrex::IntVect>& ref_ratio,
                                   amrex::Vector<amrex::Real>&          dt_level,
                                   amrex::Real                         stop_time) override;

    virtual void computeNewDt (int                                 finest_level,
                               int                                 sub_cycle,
                               amrex::Vector<int>&                  n_cycle,
                               const amrex::Vector<amrex::IntVect>& ref_ratio,
                               amrex::Vector<amrex::Real>&          dt_min,
                               amrex::Vector<amrex::Real>&          dt_level,
                               amrex::Real                         stop_time,
                               int                                 post_regrid_flag) override;

    virtual void post_regrid (int lbase, int new_finest) override;

    // Do work after timestep().
    virtual void post_timestep (int iteration) override;

    // After a full time step
    virtual void postCoarseTimeStep (amrex::Real time) override;

    // Do work after init().
    virtual void post_init (amrex::Real stop_time) override;

    virtual void post_restart () override;

    // Error estimation for regridding.
    virtual void errorEst (amrex::TagBoxArray& tb,
                           int                 clearval,
                           int                 tagval,
                           amrex::Real         time,
                           int                 n_error_buf = 0,
                           int                 ngrow = 0) override;

    virtual int WorkEstType () override { return Cost_Type; }

    // Define data descriptors.
    static void variableSetUp ();

    // Cleanup data descriptors at end of run.
    static void variableCleanUp ();

    static int numGrow() { return NUM_GROW; }

    static int numState () { return NUM_STATE; }

protected:

    static void read_params ();

    const amrex::MultiFab& volFrac () const { return *volfrac; }

    CNS& getLevel (int lev) { return dynamic_cast<CNS&>(parent->getLevel(lev)); }

    void avgDown ();

    void buildMetrics ();

    // Compute initial time step.
    amrex::Real initialTimeStep ();

    void compute_dSdt (const amrex::MultiFab& S, amrex::MultiFab& dSdt, amrex::Real dt,
                       amrex::EBFluxRegister* fr_as_crse, amrex::EBFluxRegister* fr_as_fine);

    void printTotal () const;

    const amrex::MultiFab* volfrac;
    const amrex::MultiCutFab* bndrycent;
    std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM> areafrac;
    std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM> facecent;

    amrex::iMultiFab level_mask;

    amrex::EBFluxRegister flux_reg;

    static constexpr int NUM_GROW = 5;

    enum StateDataType {
        State_Type = 0,
        Cost_Type
    };
    static int num_state_data_types;

    static amrex::BCRec phys_bc;

    // Parameters
    static int verbose;
    static amrex::IntVect hydro_tile_size;
    static amrex::Real cfl;

    static bool do_visc;
    static bool use_const_visc;

    static int plm_iorder;
    static amrex::Real plm_theta;

    static int eb_weights_type;
    static int do_reredistribution;

    static int do_reflux;

    static int refine_cutcells;

    static int refine_max_dengrad_lev;
    static amrex::Real refine_dengrad;

    static amrex::Vector<amrex::RealBox> refine_boxes;
    static amrex::RealBox* dp_refine_boxes;

    static amrex::Real gravity;

public:

    amrex::Real estTimeStep ();

    void computeTemp (amrex::MultiFab& State, int ng);

    void compute_dSdt_box (const amrex::Box& bx,
                           amrex::Array4<amrex::Real const>& Sfab,
                           amrex::Array4<amrex::Real      >& dSdtfab,
                           const std::array<FArrayBox*, AMREX_SPACEDIM>& flux);

    void compute_dSdt_box_eb (const amrex::Box& bx,
                              amrex::Array4<amrex::Real       const> const& Sfab,
                              amrex::Array4<amrex::Real            > const& dSdtfab,
                              std::array<FArrayBox*, AMREX_SPACEDIM> const& flux,
                              amrex::Array4<amrex::EBCellFlag const> const& flag,
                              amrex::Array4<amrex::Real       const> const& vfrac,
                              AMREX_D_DECL(
                              amrex::Array4<amrex::Real       const> const& apx,
                              amrex::Array4<amrex::Real       const> const& apy,
                              amrex::Array4<amrex::Real       const> const& apz),
                              AMREX_D_DECL(
                              amrex::Array4<amrex::Real       const> const& fcx,
                              amrex::Array4<amrex::Real       const> const& fcy,
                              amrex::Array4<amrex::Real       const> const& fcz),
                              amrex::Array4<amrex::Real       const> const& bcent,
                              int as_crse,
                              amrex::Array4<amrex::Real            > const& drho_as_crse,
                              amrex::Array4<int               const> const& rrflag_as_crse,
                              int as_fine,
                              amrex::Array4<amrex::Real            > const& dm_as_fine,
                              amrex::Array4<int               const> const& lev_mask,
                              amrex::Real dt);

    void cns_flux_redistribute (const amrex::Box& bx,
                                amrex::Array4<amrex::Real            > const& dqdt,
                                amrex::Array4<amrex::Real            > const& divc,
                                amrex::Array4<amrex::Real            > const& optmp,
                                amrex::Array4<amrex::Real            > const& del_m,
                                amrex::Array4<amrex::Real       const> const& redistwgt,
                                amrex::Array4<amrex::Real       const> const& vfrac,
                                amrex::Array4<amrex::EBCellFlag const> const& flag,
                                int as_crse,
                                amrex::Array4<amrex::Real            > const& drho_as_crse,
                                amrex::Array4<int               const> const& rrflag_as_crse,
                                int as_fine,
                                amrex::Array4<amrex::Real            > const& dm_as_fine,
                                amrex::Array4<int               const> const& lev_mask,
                                amrex::Real dt);

    static Parm* h_parm;
    static Parm* d_parm;
    static ProbParm* h_prob_parm;
    static ProbParm* d_prob_parm;
};

void cns_bcfill (amrex::Box const& bx, amrex::FArrayBox& data,
                 const int dcomp, const int numcomp,
                 amrex::Geometry const& geom, const amrex::Real time,
                 const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                 const int scomp);

#endif
