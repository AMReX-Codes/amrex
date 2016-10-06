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

private:

    void DepositCharge ();
    void DepositCurrent (Real dt);
    void EvolveE (Real dt);
    void EvolveB (Real dt);

    void ReadParameters ();

    //
    // data
    //
    Array<Geometry> geom_arr;
    Array<BoxArray> ba_arr;
    Array<DistributionMapping> dmap_arr;
    // Particle container
    std::unique_ptr<MyParticleContainer> mypc;
    // Assuming there is only a single level
    MultiFab charge;
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
    int nx, ny, nz;
    int max_step;
    //
    int charge_deposition_order = 1;
};

#endif
