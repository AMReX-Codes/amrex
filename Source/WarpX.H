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

#include <ParticleContainer.H>
#include <WarpXCrseFineBndry.H>

/*! A test class */

class WarpX
    : public amrex::AmrCore
{
public:

    static WarpX& GetInstance ();
    static void ResetInstance ();

    WarpX ();
    ~WarpX ();

    int Verbose () const { return verbose; }
    int RegridInt () const { return regrid_int; }

    void InitData ();

    void Evolve (int numsteps = -1);

    MultiParticleContainer& GetPartContainer () { return *mypc; }

    static void shiftMF(amrex::MultiFab& mf, const amrex::Geometry& geom, int num_shift, int dir);

    static void GotoNextLine (std::istream& is);

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

    const amrex::MultiFab& getcurrent (int lev, int direction) {return *current[lev][direction];}
    const amrex::MultiFab& getEfield (int lev, int direction) {return *Efield[lev][direction];}
    const amrex::MultiFab& getBfield (int lev, int direction) {return *Bfield[lev][direction];}

    void ComputeDt ();
    void MoveWindow (bool move_j);

    const amrex::Array<amrex::Real>& getPMLSigma(int dir) {return pml_sigma[dir];};
    const amrex::Array<amrex::Real>& getPMLSigmaStar(int dir) {return pml_sigma_star[dir];};
    void ComputePMLFactors (int lev, amrex::Real dt);
    void ComputePMLFactorsE (int lev, amrex::Real dt);
    void ComputePMLFactorsB (int lev, amrex::Real dt);
    void ExchangeWithPML (amrex::MultiFab& regmf, amrex::MultiFab& pmlmf,
                          const amrex::Geometry& gm);

    void EvolveE (         amrex::Real dt);
    void EvolveE (int lev, amrex::Real dt);
    void EvolveB (         amrex::Real dt);
    void EvolveB (int lev, amrex::Real dt);

    void PushParticlesandDepose (int lev, amrex::Real cur_time);

    // Fill boundary cells including coarse/fine boundaries
    void FillBoundaryB ();
    void FillBoundaryE ();
    void FillBoundaryE (int lev, bool force);
    void FillBoundaryB (int lev, bool force);

    void AverageDownB ();
    void AverageDownE ();
    void AverageDownCurrent ();
    void SyncCurrent ();

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
			      const amrex::DistributionMapping& dm) final
        { amrex::Abort("RemakeLevel: To be implemented"); }

    //! Delete level data.  Called by AmrCore::regrid.
    virtual void ClearLevel (int lev) final;

private:

    // Singleton is used when the code is run from python
    static WarpX* m_instance;

    void ReadParameters ();

    // *************************************************************************
    void RegridBaseLevel ();
    bool LoadBalanceBaseLevel ();
    void splitBoxes (amrex::BoxArray& ba, amrex::Array<long>& newcost, const amrex::Array<long>& cost_in, int heavy_grid_size);
    amrex::Real getEfficiency(const amrex::DistributionMapping& dm, const amrex::Array<long>& cost);
    amrex::DistributionMapping getCostCountDM (const amrex::Array<long>& cost, const amrex::BoxArray& ba);
    bool okToRegrid(int step);
    // *************************************************************************

    void InitFromScratch ();

    void AllocLevelData (int lev, const amrex::BoxArray& new_grids,
                         const amrex::DistributionMapping& new_dmap);
    void InitLevelData (int lev, amrex::Real time);

    void InitFromCheckpoint ();
    void PostRestart ();

    void InitOpenbc ();

    void InjectPlasma(int num_shift, int dir);

    void WriteWarpXHeader(const std::string& name) const;
    void WriteJobInfo (const std::string& dir) const;
    static void PackPlotDataPtrs (amrex::Array<const amrex::MultiFab*>& pmf,
				  const std::array<std::unique_ptr<amrex::MultiFab>,3>& data);

    void InitPML ();

    static void ComputeDivB (amrex::MultiFab& divB, int dcomp,
                             const amrex::Array<const amrex::MultiFab*>& B,
                             const amrex::Real* dx);

    void AverageDownEorCurrent (amrex::Array<std::array<std::unique_ptr<amrex::MultiFab>,3> >& data);

    amrex::Array<int> istep;      // which step?
    amrex::Array<int> nsubsteps;  // how many substeps on each level?

    amrex::Array<amrex::Real> t_new;
    amrex::Array<amrex::Real> t_old;
    amrex::Array<amrex::Real> dt;

    // Particle container
    std::unique_ptr<MultiParticleContainer> mypc;

    // Fields: First array for level, second for direction
    amrex::Array<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current;
    amrex::Array<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield;
    amrex::Array<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield;

    amrex::Array<std::unique_ptr<WarpXCrseFineBndry> > cfbndry;

    // PML  xxxxx for now let's only consider single level
    int do_pml = 0;
    int pml_ncell = 10;
    amrex::BoxArray pml_ba;
    amrex::DistributionMapping pml_dm;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_E;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_B;
    //
    struct PMLSigma : amrex::Array<amrex::Real>
    {
        int lo() const { return m_lo; }
        int hi() const { return m_hi; }
        int m_lo, m_hi;
    };
    std::array<PMLSigma,BL_SPACEDIM> pml_sigma, pml_sigma_star; // sigma/epsilon, sigma_star/mu
    std::array<PMLSigma,BL_SPACEDIM> pml_sigma_fac1, pml_sigma_fac2, 
                                     pml_sigma_star_fac1, pml_sigma_star_fac2;

    // Moving window parameters
    int do_moving_window  = 0;
    int moving_window_dir = -1;
    amrex::Real moving_window_x  = std::numeric_limits<amrex::Real>::max();
    amrex::Real moving_window_v  = std::numeric_limits<amrex::Real>::max();

    // Plasma injection parameters
    int do_plasma_injection  = 0;
    int num_injected_species = -1;
    amrex::Array<int> injected_plasma_species;

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
    bool plot_part_per_grid = true;
    bool plot_part_per_proc = true;
    bool plot_proc_number   = true;
    bool plot_divb          = false;
    bool plot_raw_fields = false;

    bool is_synchronized = true;
};

#endif
