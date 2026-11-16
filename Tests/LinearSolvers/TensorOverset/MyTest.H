#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MLMG.H>

class MyTest
{
public:

    MyTest ();

    void solve ();
    void writePlotfile ();
    void initData ();

private:

    void readParameters ();
    void initGrids ();

    int n_cell = 128;
    int max_grid_size = 64;

    std::string plot_file_name{"plot"};

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 2;
    int max_coarsening_level = 30;

    amrex::Geometry geom;
    amrex::BoxArray grids;
    amrex::DistributionMapping dmap;

    amrex::MultiFab solution;
    amrex::MultiFab exact;
    amrex::MultiFab rhs;
    amrex::MultiFab eta;

    int do_overset = 1;
    int mixed_precision = 0;
    amrex::Real max_error = -1.0; // if positive, abort when the max-norm error exceeds it
    int do_plots = 1;
    amrex::iMultiFab oversetmask;
};

#endif
