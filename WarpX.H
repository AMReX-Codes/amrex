#ifndef WARPX_H_
#define WARPX_H_

#include <iostream>
#include <memory>

#include <Array.H>
#include <BLProfiler.H>
#include <Geometry.H>
#include <MultiFab.H>

#include <ParticleContainer.H>

class WarpX
{
public:
    WarpX ();
    ~WarpX ();

    void InitData ();

    void Evolve ();

    int Verbose () const { return verbose; }

    static void FillBoundary (MultiFab& mf, const Geometry& geom, const IntVect& nodalflag);

private:

    void EvolveE (Real dt);
    void EvolveB (Real dt);

    void ReadParameters ();

    void WritePlotFile (int istep, Real t) const;

    //
    // data
    //
    Array<Geometry> geom_arr;
    Array<BoxArray> ba_arr;
    Array<DistributionMapping> dmap_arr;
    //
    // Particle container
    std::unique_ptr<MyParticleContainer> mypc;
    //
    // Assuming there is only a single level
    Array< std::unique_ptr<MultiFab> > current;
    Array< std::unique_ptr<MultiFab> > Efield;
    Array< std::unique_ptr<MultiFab> > Bfield;

    //
    // parameters
    //
    int verbose = 1;
    //
    int max_level = 0;
    int max_grid_size = 32;
    int ref_ratio = 2;
    //
    Real cfl = 1;
    int nx, ny, nz;
    int max_step;
    //
    int plot_int = -1;
};

#endif
