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

    static WarpX& GetInstance ();
    static void ResetInstance ();

    WarpX ();
    ~WarpX ();

    int Verbose () const { return verbose; }
    int RegridInt () const { return regrid_int; }

    void InitData ();

    void Evolve (int numsteps = -1);

    MultiParticleContainer& GetPartContainer () { return *mypc; }

    void fillSlice(Real z_coord) const;

    void sampleAtPoints(const Array<Real>& x,
			const Array<Real>& y,
			const Array<Real>& z,
			Array<Array<Real> >& result) const;

    static void FillBoundary (MultiFab& mf, const Geometry& geom, const IntVect& nodalflag);

    static void shiftMF(MultiFab& mf, const Geometry& geom, int num_shift, 
			int dir, const IntVect& nodalflag);

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

protected:
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/) override {}

    Box getIndexBox(const RealBox& real_box) const;

private:

    // Singleton is used when the code is run from python
    static WarpX* m_instance;

    void EvolveE (int lev, Real dt);
    void EvolveB (int lev, Real dt);

    void ReadParameters ();

    void MakeNewLevel (int lev, const BoxArray& new_grids, const DistributionMapping& new_dmap);

    void RegridBaseLevel ();
    bool LoadBalanceBaseLevel ();
    void splitBoxes (BoxArray& ba, Array<long>& newcost, const Array<long>& cost_in, int heavy_grid_size);

    Real getEfficiency(const DistributionMapping& dm, const Array<long>& cost);
    DistributionMapping getCostCountDM (const Array<long>& cost, const BoxArray& ba);

    void InitFromScratch ();
    void InitFromCheckpoint ();
    void InitLevelData (int lev);
    void PostRestart ();

    void InitOpenbc ();

    void ComputeDt ();

    void MoveWindow ();
    void InjectPlasma(int num_shift, int dir);

    void WriteWarpXHeader(const std::string& name) const;
    void WriteCheckPointFile() const;
    void WritePlotFile () const;
    void WriteJobInfo (const std::string& dir) const;
    static void PackPlotDataPtrs (std::vector<MultiFab*>& pmf,
				  const Array<std::unique_ptr<MultiFab> >& data);
    static void Copy (MultiFab& dstmf, int dcomp, int ncomp, const MultiFab& srcmf, int scomp);

    bool okToRegrid(int step);

    Array<int> istep;      // which step?
    Array<int> nsubsteps;  // how many substeps on each level?

    Array<Real> t_new;
    Array<Real> t_old;
    Array<Real> dt;

    // Particle container
    std::unique_ptr<MultiParticleContainer> mypc;

    // Fields: First array for level, second for direction
    Array<Array< std::unique_ptr<MultiFab> > > current;
    Array<Array< std::unique_ptr<MultiFab> > > Efield;
    Array<Array< std::unique_ptr<MultiFab> > > Bfield;

    // Moving window parameters
    int do_moving_window  = 0;
    int moving_window_dir = -1;
    Real moving_window_x  = std::numeric_limits<Real>::max();
    Real moving_window_v  = std::numeric_limits<Real>::max();

    // Plasma injection parameters
    int do_plasma_injection  = 0;
    int num_injected_species = -1;
    Array<int> injected_plasma_ppc;
    Array<int> injected_plasma_species;
    Array<Real> injected_plasma_density;

    // Other runtime parameters
    int verbose = 1;

    int max_step   = std::numeric_limits<int>::max();
    Real stop_time = std::numeric_limits<Real>::max();

    int regrid_int = -1;

    Real cfl = 0.7;

    std::string restart_chkfile;

    std::string check_file {"chk"};
    std::string plot_file {"plt"};
    int check_int = -1;
    int plot_int = -1;

    static IntVect Bx_nodal_flag;
    static IntVect By_nodal_flag;
    static IntVect Bz_nodal_flag;

    static IntVect Ex_nodal_flag;
    static IntVect Ey_nodal_flag;
    static IntVect Ez_nodal_flag;

    static IntVect jx_nodal_flag;
    static IntVect jy_nodal_flag;
    static IntVect jz_nodal_flag;
};

#endif
