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
    void checkOversetCells () const; // public for CUDA device lambdas

private:

    void readParameters ();
    void initGrids ();

    int max_level = 0;
    int n_cell = 128;
    int max_grid_size = 64;

    std::string plot_file_name{"plot"};

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 2;
    int max_coarsening_level = 30;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    amrex::Vector<amrex::MultiFab> phi;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> exact_phi;
    amrex::Vector<amrex::MultiFab> acoef;
    amrex::Vector<amrex::MultiFab> bcoef;

    amrex::Real ascalar = 1.e-3;
    amrex::Real bscalar = 1.0;

    int do_overset = 1;
    amrex::Vector<amrex::iMultiFab> oversetmask;

    int use_hypre = 0;
};

#endif
