#ifndef WARPX_H_
#define WARPX_H_

#include <iostream>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_BLProfiler.H>

#include <ParticleContainer.H>

class WarpX
    : public amrex::AmrCore
{
public:

    static WarpX& GetInstance ();
    static void ResetInstance ();

    WarpX ();
    ~WarpX ();

    int Verbose () const { return verbose; }

    void InitData ();

    void Evolve (int numsteps = -1);

    static void FillBoundary (amrex::MultiFab& mf, const amrex::Geometry& geom,
			      const amrex::IntVect& nodalflag);

    // Algorithms
    static long current_deposition_algo;
    static long charge_deposition_algo;
    static long field_gathering_algo;
    static long particle_pusher_algo;
    
    // Interpolation order
    static long nox;
    static long noy;
    static long noz;

protected:
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/) override {}

private:

    // Singleton is used when the code is run from python
    static WarpX* m_instance;

    void EvolveE (int lev, amrex::Real dt);
    void EvolveB (int lev, amrex::Real dt);

    void ReadParameters ();

    void MakeNewLevel (int lev, amrex::Real time,
		       const amrex::BoxArray& new_grids, const amrex::DistributionMapping& new_dmap);

    void InitFromScratch ();
    void InitFromCheckpoint () {}
    void InitLevelData (int lev);
    
    void InitOpenbc ();
    
    void ComputeDt ();

    void WritePlotFile () const;
    static void PackPlotDataPtrs (amrex::Array<const amrex::MultiFab*>& pmf,
				  const amrex::Array<std::unique_ptr<amrex::MultiFab> >& data);
    static void Copy (amrex::MultiFab& dstmf, int dcomp, int ncomp, const amrex::MultiFab& srcmf, int scomp);

    amrex::Array<int> istep;      // which step?
    amrex::Array<int> nsubsteps;  // how many substeps on each level?

    amrex::Array<amrex::Real> t_new;
    amrex::Array<amrex::Real> t_old;
    amrex::Array<amrex::Real> dt;

    // Particle container
    std::unique_ptr<MyParticleContainer> mypc;
    // Fields: First array for level, second for direction
    amrex::Array<amrex::Array< std::unique_ptr<amrex::MultiFab> > > current;
    amrex::Array<amrex::Array< std::unique_ptr<amrex::MultiFab> > > Efield;
    amrex::Array<amrex::Array< std::unique_ptr<amrex::MultiFab> > > Bfield;

    // runtime parameters

    int verbose = 1;

    int max_step   = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    amrex::Real cfl = 0.7;

    std::string restart_chkfile;

    std::string check_file {"chk"};
    std::string plot_file {"plt"};
    int check_int = -1;
    int plot_int = -1;

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
