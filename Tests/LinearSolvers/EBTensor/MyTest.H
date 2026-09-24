#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MLMG.H>

#ifdef AMREX_USE_HYPRE
#include <AMReX_Hypre.H>
#endif

class MyTest
{
public:

    MyTest ();

    void solve ();
    void writePlotfile () const;

    // public for cuda
    void initProb ();
    void initData ();

private:

    void readParameters ();

    int n_cell = 128;
    amrex::IntVect n_cells; // number of cells in each direction, n_cell by default
    int max_grid_size = 64;

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 1;
    int max_iter = 100;
    int max_fmg_iter = 0;
    int linop_maxorder = 2;
    bool agglomeration = true;
    bool consolidation = true;
    int max_coarsening_level = 30;
    amrex::Real max_error = -1; // abort if a max-norm error exceeds it (if > 0)

    amrex::Geometry geom;
    amrex::BoxArray grids;
    amrex::DistributionMapping dmap;
    std::unique_ptr<amrex::EBFArrayBoxFactory> factory;

    amrex::MultiFab solution;
    amrex::MultiFab exact;
    amrex::MultiFab rhs;
    amrex::MultiFab eta;
};

#endif
