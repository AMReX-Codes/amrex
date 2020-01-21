#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLMG.H>

class MyTest
{
public:

    MyTest ();

    void solve ();
    void writePlotfile ();
    void initData ();

private:

    void initializeEB ();
    void readParameters ();
    void initGrids ();

    int max_level = 0;
    int ref_ratio = 2;
    int n_cell = 128;
    int max_grid_size = 64;

    std::string plot_file_name{"plot"};

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 2;
    int max_iter = 100;
    int max_fmg_iter = 0;
    int max_coarsening_level = 30;
    bool use_hypre = false;
    int agg_grid_size = -1;
    int con_grid_size = -1;

    int cylinder_direction = 0;

    bool sigma = false;  // only for all_regular, we can set it to true

    int gpu_regtest = 0;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;
    amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory> > factory;

    amrex::Vector<amrex::MultiFab> phi;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> vel;
    amrex::Vector<amrex::MultiFab> sig;
};

#endif
