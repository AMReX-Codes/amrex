#ifndef WARPX_H_
#define WARPX_H_

#include <iostream>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>

#include <ParticleContainer.H>

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

    static void FillBoundary (amrex::MultiFab& mf, const amrex::Geometry& geom, const amrex::IntVect& nodalflag);

    static void shiftMF(amrex::MultiFab& mf, const amrex::Geometry& geom, int num_shift,
			int dir, const amrex::IntVect& nodalflag);

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

    amrex::MultiFab * getcurrent (int lev, int direction) {return current[lev][direction].get();}
    amrex::MultiFab * getEfield (int lev, int direction) {return Efield[lev][direction].get();}
    amrex::MultiFab * getBfield (int lev, int direction) {return Bfield[lev][direction].get();}

    void ComputeDt ();
    void MoveWindow ();

    void EvolveE (int lev, amrex::Real dt);
    void EvolveB (int lev, amrex::Real dt);
    void FillBoundaryE (int lev, bool force);
    void FillBoundaryB (int lev, bool force);
    void PushParticlesandDepose (int lev, amrex::Real cur_time);

    int getistep (int lev) {return istep[lev];}
    void setistep (int lev, int ii) {istep[lev] = ii;}
    amrex::Real gett_new (int lev) {return t_new[lev];}
    void sett_new (int lev, amrex::Real time) {t_new[lev] = time;}
    amrex::Real getdt (int lev) {return dt[lev];}

    int maxStep () {return max_step;}
    amrex::Real stopTime () {return stop_time;}

    int checkInt () {return check_int;}
    int plotInt () {return plot_int;}

    void WriteCheckPointFile () const;
    void WritePlotFile () const;

protected:

    //! Tagging cells for refinement
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/) override
        { amrex::Abort("ErrorEst: To be implemented"); }

    //! Make a new level from scratch using provided BoxArray and
    //! DistributionMapping.  Only used during initialization.  Called
    //! by AmrCoreInitFromScratch.
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					  const amrex::DistributionMapping& dm) override;

    //! Make a new level using provided BoxArray and
    //! DistributionMapping and fill with interpolated coarse level
    //! data.  Called by AmrCore::regrid.
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
					 const amrex::DistributionMapping& dm) override
        { amrex::Abort("MakeNewLevelFromCoarse: To be implemented"); }

    //! Remake an existing level using provided BoxArray and
    //! DistributionMapping and fill with existing fine and coarse
    //! data.  Called by AmrCore::regrid.
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
			      const amrex::DistributionMapping& dm) override
        { amrex::Abort("RemakeLevel: To be implemented"); }

    //! Delete level data.  Called by AmrCore::regrid.
    virtual void ClearLevel (int lev) override;

private:

    // Singleton is used when the code is run from python
    static WarpX* m_instance;

    void ReadParameters ();

    void RegridBaseLevel ();
    bool LoadBalanceBaseLevel ();
    void splitBoxes (amrex::BoxArray& ba, amrex::Array<long>& newcost, const amrex::Array<long>& cost_in, int heavy_grid_size);

    amrex::Real getEfficiency(const amrex::DistributionMapping& dm, const amrex::Array<long>& cost);
    amrex::DistributionMapping getCostCountDM (const amrex::Array<long>& cost, const amrex::BoxArray& ba);

    void InitFromScratch ();

    void AllocLevelData (int lev, const amrex::BoxArray& new_grids,
                         const amrex::DistributionMapping& new_dmap);
    void InitLevelData (int lev);

    void InitFromCheckpoint ();
    void PostRestart ();

    void InitOpenbc ();

    void InjectPlasma(int num_shift, int dir);

    void WriteWarpXHeader(const std::string& name) const;
    void WriteJobInfo (const std::string& dir) const;
    static void PackPlotDataPtrs (amrex::Array<const amrex::MultiFab*>& pmf,
				  const amrex::Array<std::unique_ptr<amrex::MultiFab> >& data);
    static void Copy (amrex::MultiFab& dstmf, int dcomp, int ncomp, const amrex::MultiFab& srcmf, int scomp);

    bool okToRegrid(int step);

    amrex::Array<int> istep;      // which step?
    amrex::Array<int> nsubsteps;  // how many substeps on each level?

    amrex::Array<amrex::Real> t_new;
    amrex::Array<amrex::Real> t_old;
    amrex::Array<amrex::Real> dt;

    // Particle container
    std::unique_ptr<MultiParticleContainer> mypc;
    // Fields: First array for level, second for direction
    amrex::Array<amrex::Array< std::unique_ptr<amrex::MultiFab> > > current;
    amrex::Array<amrex::Array< std::unique_ptr<amrex::MultiFab> > > Efield;
    amrex::Array<amrex::Array< std::unique_ptr<amrex::MultiFab> > > Bfield;

    // Moving window parameters
    int do_moving_window  = 0;
    int moving_window_dir = -1;
    amrex::Real moving_window_x  = std::numeric_limits<amrex::Real>::max();
    amrex::Real moving_window_v  = std::numeric_limits<amrex::Real>::max();

    // Plasma injection parameters
    int do_plasma_injection  = 0;
    int num_injected_species = -1;
    amrex::Array<int> injected_plasma_ppc;
    amrex::Array<int> injected_plasma_species;
    amrex::Array<amrex::Real> injected_plasma_density;

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

    bool is_synchronized;

    static amrex::IntVect Bx_nodal_flag;
    static amrex::IntVect By_nodal_flag;
    static amrex::IntVect Bz_nodal_flag;

    static amrex::IntVect Ex_nodal_flag;
    static amrex::IntVect Ey_nodal_flag;
    static amrex::IntVect Ez_nodal_flag;

    static amrex::IntVect jx_nodal_flag;
    static amrex::IntVect jy_nodal_flag;
    static amrex::IntVect jz_nodal_flag;
};

#endif
