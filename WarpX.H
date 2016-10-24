#ifndef WARPX_H_
#define WARPX_H_

#include <iostream>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AmrCore.H>
#include <BLProfiler.H>

#include <ParticleContainer.H>

class WarpX
    : public AmrCore
{
public:
    WarpX ();
    ~WarpX ();

    int Verbose () const { return verbose; }

    void InitData ();

    void Evolve ();

    static void FillBoundary (MultiFab& mf, const Geometry& geom, const IntVect& nodalflag);

protected:
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/) override {}

private:

    void EvolveE (int lev, Real dt);
    void EvolveB (int lev, Real dt);

    void ReadParameters ();

    void MakeNewLevel (int lev, Real time,
		       const BoxArray& new_grids, const DistributionMapping& new_dmap);

    void InitFromScratch ();
    void InitFromCheckpoint () {}
    void InitLevelData (int lev);
    
    void ComputeDt ();

    void WritePlotFile () const;

    Array<int> istep;      // which step?
    Array<int> nsubsteps;  // how many substeps on each level?

    Array<Real> t_new;
    Array<Real> t_old;
    Array<Real> dt;

    // Particle container
    std::unique_ptr<MyParticleContainer> mypc;
    // Fields: First array for level, second for direction
    Array<Array< std::unique_ptr<MultiFab> > > current;
    Array<Array< std::unique_ptr<MultiFab> > > Efield;
    Array<Array< std::unique_ptr<MultiFab> > > Bfield;

    // runtime parameters

    int verbose = 1;

    int max_step   = std::numeric_limits<int>::max();
    Real stop_time = std::numeric_limits<Real>::max();

    Real cfl = 0.7;

    std::string restart_chkfile;

    std::string check_file {"chk"};
    std::string plot_file {"plt"};
    int check_int = -1;
    int plot_int = -1;
};

#endif
