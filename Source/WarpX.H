#ifndef WARPX_H_
#define WARPX_H_

#include <iostream>
#include <memory>
#include <array>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>
#include <AMReX_RealVect.H>
#include <AMReX_VisMF.H>

#include <ParticleContainer.H>
#include <WarpXPML.H>

class NoOpPhysBC
    : public amrex::PhysBCFunctBase
{
public:
    NoOpPhysBC () {}
    virtual ~NoOpPhysBC () {}
    virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override { }
    using amrex::PhysBCFunctBase::FillBoundary;
};

enum struct DtType : int
{
    Full = 0,
    FirstHalf,
    SecondHalf
};

class WarpX
    : public amrex::AmrCore
{
public:

    friend class PML;

    static WarpX& GetInstance ();
    static void ResetInstance ();

    WarpX ();
    ~WarpX ();

    int Verbose () const { return verbose; }

    void InitData ();

    void Evolve (int numsteps = -1);

    MultiParticleContainer& GetPartContainer () { return *mypc; }

    static void shiftMF(amrex::MultiFab& mf, const amrex::Geometry& geom, int num_shift, int dir);

    static void GotoNextLine (std::istream& is);

    // External fields
    static amrex::Vector<amrex::Real> B_external;

    // Algorithms
    static long current_deposition_algo;
    static long charge_deposition_algo;
    static long field_gathering_algo;
    static long particle_pusher_algo;

    // Interpolation order
    static long nox;
    static long noy;
    static long noz;

    static bool use_laser;
    static bool use_filter;
    static bool serialize_ics;

    // Boosted frame parameters
    static amrex::Real gamma_boost;
    static amrex::Real beta_boost;
    static amrex::Vector<amrex::Real> boost_direction;

    const amrex::MultiFab& getcurrent (int lev, int direction) {return *current_fp[lev][direction];}
    const amrex::MultiFab& getEfield (int lev, int direction) {return *Efield_fp[lev][direction];}
    const amrex::MultiFab& getBfield (int lev, int direction) {return *Bfield_fp[lev][direction];}

    static amrex::MultiFab* getCosts (int lev) {
        if (m_instance) {
            return m_instance->costs[lev].get();
        } else {
            return nullptr;
        }
    }

    void ComputeDt ();
    void MoveWindow (bool move_j);

    void EvolveE (         amrex::Real dt, DtType typ);
    void EvolveE (int lev, amrex::Real dt, DtType typ);
    void EvolveB (         amrex::Real dt, DtType typ);
    void EvolveB (int lev, amrex::Real dt, DtType typ);
    void EvolveF (         amrex::Real dt, DtType typ);
    void EvolveF (int lev, amrex::Real dt, DtType typ);

    void PushParticlesandDepose (int lev, amrex::Real cur_time);
    void PushParticlesandDepose (         amrex::Real cur_time);

    // This function does aux(lev) = fp(lev) + I(aux(lev-1)-cp(lev)).
    // Caller must make sure fp and cp have ghost cells filled.
    void UpdateAuxilaryData ();

    // Fill boundary cells including coarse/fine boundaries
    void FillBoundaryB ();
    void FillBoundaryE ();
    void FillBoundaryE (int lev);
    void FillBoundaryB (int lev);

    void SyncCurrent ();
    void SyncRho ();

    int getistep (int lev) const {return istep[lev];}
    void setistep (int lev, int ii) {istep[lev] = ii;}
    amrex::Real gett_new (int lev) const {return t_new[lev];}
    void sett_new (int lev, amrex::Real time) {t_new[lev] = time;}
    amrex::Real getdt (int lev) const {return dt[lev];}

    int maxStep () const {return max_step;}
    amrex::Real stopTime () const {return stop_time;}

    int checkInt () const {return check_int;}
    int plotInt () const {return plot_int;}

    void WriteCheckPointFile () const;
    void WritePlotFile () const;

    void WritePlotFileES(const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                         const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                         const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E);

    static std::array<amrex::Real,3> CellSize (int lev);
    static std::array<amrex::Real,3> LowerCorner (const amrex::Box& bx, int lev);
    static std::array<amrex::Real,3> UpperCorner (const amrex::Box& bx, int lev);

    static amrex::IntVect Bx_nodal_flag;
    static amrex::IntVect By_nodal_flag;
    static amrex::IntVect Bz_nodal_flag;

    static amrex::IntVect Ex_nodal_flag;
    static amrex::IntVect Ey_nodal_flag;
    static amrex::IntVect Ez_nodal_flag;

    static amrex::IntVect jx_nodal_flag;
    static amrex::IntVect jy_nodal_flag;
    static amrex::IntVect jz_nodal_flag;

protected:

    //! Tagging cells for refinement
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/) final;

    //! Make a new level from scratch using provided BoxArray and
    //! DistributionMapping.  Only used during initialization.  Called
    //! by AmrCoreInitFromScratch.
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					  const amrex::DistributionMapping& dm) final;

    //! Make a new level using provided BoxArray and
    //! DistributionMapping and fill with interpolated coarse level
    //! data.  Called by AmrCore::regrid.
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
					 const amrex::DistributionMapping& dm) final
        { amrex::Abort("MakeNewLevelFromCoarse: To be implemented"); }

    //! Remake an existing level using provided BoxArray and
    //! DistributionMapping and fill with existing fine and coarse
    //! data.  Called by AmrCore::regrid.
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
			      const amrex::DistributionMapping& dm) final;

    //! Delete level data.  Called by AmrCore::regrid.
    virtual void ClearLevel (int lev) final;

private:

    // Singleton is used when the code is run from python
    static WarpX* m_instance;

    ///
    /// Advance the simulation by numsteps steps, electromagnetic case.
    ///
    void EvolveEM(int numsteps);

    ///
    /// Advance the simulation by numsteps steps, electrostatic case.
    ///
    void EvolveES(int numsteps);

    ///
    /// Compute the gravitational potential from rho by solving Poisson's equation.
    /// Both rho and phi are assumed to be node-centered. This method is only used
    /// in electrostatic mode.
    ///
    void computePhi(const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                          amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi) const;

    ///
    /// Compute the electric field in each direction by computing the gradient
    /// the potential phi using 2nd order centered differences. Both rho and phi
    /// are assumed to be node-centered. This method is only used in electrostatic mode.
    ///
    void computeE(amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                  const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi) const;

    void ReadParameters ();

    void InitFromScratch ();

    void AllocLevelData (int lev, const amrex::BoxArray& new_grids,
                         const amrex::DistributionMapping& new_dmap);
    void InitLevelData (int lev, amrex::Real time);

    void InitFromCheckpoint ();
    void PostRestart ();

    void InitOpenbc ();

    void InitPML ();
    void ComputePMLFactors ();

    void InjectPlasma(int num_shift, int dir);

    void WriteWarpXHeader(const std::string& name) const;
    void WriteJobInfo (const std::string& dir) const;
    static void PackPlotDataPtrs (amrex::Vector<const amrex::MultiFab*>& pmf,
				  const std::array<std::unique_ptr<amrex::MultiFab>,3>& data);

    static void ComputeDivB (amrex::MultiFab& divB, int dcomp,
                             const std::array<const amrex::MultiFab*, 3>& B,
                             const std::array<amrex::Real,3>& dx);

    static void ComputeDivE (amrex::MultiFab& divE, int dcomp,
                             const std::array<const amrex::MultiFab*, 3>& B,
                             const std::array<amrex::Real,3>& dx);

    void SyncCurrent (const std::array<const amrex::MultiFab*,3>& fine,
                      const std::array<      amrex::MultiFab*,3>& crse,
                      int ref_ratio);

    void SyncRho (const amrex::MultiFab& fine, amrex::MultiFab& crse, int ref_ratio);

    void ExchangeWithPmlB (int lev);
    void ExchangeWithPmlE (int lev);

    //
    // This stuff is needed by the nodal multigrid solver when running in
    // electrostatic mode.
    //
    void zeroOutBoundary(amrex::MultiFab& input_data, amrex::MultiFab& bndry_data,
                         const amrex::FabArray<amrex::BaseFab<int> >& mask) const;

    void sumFineToCrseNodal(const amrex::MultiFab& fine, amrex::MultiFab& crse,
                            const amrex::Geometry& cgeom, const amrex::IntVect& ratio);

    void fixRHSForSolve(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rhs,
                        const amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks) const ;

    void getLevelMasks(amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks,
                       const int nnodes = 1);

    void LoadBalance ();

    // used to zero out fine level data on points shared with the coarse grid
    // in electrostatic mode
    amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > > masks;

    // used to gather the field from the coarse level in electrostatic mode.
    amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > > gather_masks;

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // Particle container
    std::unique_ptr<MultiParticleContainer> mypc;

    //
    // Fields: First array for level, second for direction
    //

    // Full solution
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_aux;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_aux;

    // Fine patch
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > F_fp;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > rho_fp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_fp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_fp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_fp;

    // Coarse patch
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > F_cp;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > rho_cp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_cp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_cp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_cp;

    // div E cleaning
    int do_dive_cleaning = 0;

    // PML
    int do_pml = 1;
    int pml_ncell = 10;
    int pml_delta = 10;
    std::string pml_type = "pml";   // or apml
    amrex::Vector<std::unique_ptr<PML> > pml;

    // Moving window parameters
    int do_moving_window  = 0;
    int moving_window_dir = -1;
    amrex::Real moving_window_x  = std::numeric_limits<amrex::Real>::max();
    amrex::Real moving_window_v  = std::numeric_limits<amrex::Real>::max();

    // Plasma injection parameters
    int do_plasma_injection  = 0;
    int num_injected_species = -1;
    amrex::Vector<int> injected_plasma_species;

    int do_electrostatic = 0;
    int n_buffer = 4;
    amrex::Real const_dt = 0.5e-11;

    int load_balance_int = -1;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > costs;

    // Other runtime parameters
    int verbose = 1;

    int max_step   = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    int regrid_int = -1;

    amrex::Real cfl = 0.7;

    std::string restart_chkfile;

    std::string check_file {"chk"};
    std::string plot_file {"plt"};
    int check_int = -1;
    int plot_int = -1;

    bool plot_part_per_cell = true;
    bool plot_part_per_grid = false;
    bool plot_part_per_proc = false;
    bool plot_proc_number   = false;
    bool plot_dive          = false;
    bool plot_divb          = true;
    bool plot_finepatch     = false;
    bool plot_crsepatch     = false;
    bool plot_raw_fields    = false;
    bool plot_raw_fields_guards = false;

    amrex::VisMF::Header::Version checkpoint_headerversion = amrex::VisMF::Header::NoFabHeader_v1;
    amrex::VisMF::Header::Version plot_headerversion       = amrex::VisMF::Header::Version_v1;
    bool use_single_read = true;
    bool use_single_write = true;
    int mffile_nstreams = 4;
    int field_io_nfiles = 1024;
    int particle_io_nfiles = 1024;

    amrex::RealVect fine_tag_lo;
    amrex::RealVect fine_tag_hi;

    bool is_synchronized = true;
};

#endif
