#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MLMG.H>

class MyTest
{
public:

    MyTest ();

    void solve ();

private:

    void readParameters ();
    void initData ();

    int max_level = 0;
    int ref_ratio = 2;
    int n_cell = 128;
    int max_grid_size = 64;
    
    // For MLMG solver
    int verbose = 2;
    int cg_verbose = 2;
    int max_iter = 100;
    int max_fmg_iter = 0;
    bool use_hypre = false;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    amrex::Vector<amrex::MultiFab> phi;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> vel;
};

#endif
