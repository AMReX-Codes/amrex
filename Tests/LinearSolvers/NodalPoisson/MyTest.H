#ifndef MY_TEST_H_
#define MY_TEST_H_
#ifdef AMREX_USE_HYPRE
#include <AMReX_Hypre.H>
#endif

#include <AMReX_MLMG.H>

// Solve `div dot (sigma grad phi) = rhs` using the nodal solver

class MyTest
{
public:

    MyTest ();

    void solve ();
    void compute_norms () const;
    void writePlotfile () const;

    void initData ();

    int getNumTrials() { return num_trials; }
    int getDoPlots() { return do_plots; }

private:

    void readParameters ();

    int max_level = 1;
    int ref_ratio = 2;
    int n_cell = 128;
    int max_grid_size = 64;

    bool composite_solve = true;

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 0;
    int max_iter = 100;
    int max_fmg_iter = 0;
    amrex::Real reltol = 1.e-11;
    int gpu_regtest = 0;


    bool agglomeration = false;
    bool consolidation = false;
    bool semicoarsening = false;
    int max_coarsening_level = 30;
    int max_semicoarsening_level = 0;
    //int smooth_num_sweeps = 4;

    bool use_hypre = false;
    bool do_plots = true;
    int num_trials = 1;

#ifdef AMREX_USE_HYPRE
    int hypre_interface_i = 1;  // 1. structed, 2. semi-structed, 3. ij
    amrex::Hypre::Interface hypre_interface = amrex::Hypre::Interface::structed;
#endif

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    amrex::Vector<amrex::MultiFab> solution;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> exact_solution;
    amrex::Vector<amrex::MultiFab> sigma;
};

#endif
