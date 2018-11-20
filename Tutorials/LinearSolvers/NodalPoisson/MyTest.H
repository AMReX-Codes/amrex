#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MLMG.H>

// Solve `div dot (sigma grad phi) = rhs` using the nodal solver

class MyTest
{
public:

    MyTest ();

    void solve ();

    void compute_norms () const;

private:

    void readParameters ();
    void initData ();

    int n_cell = 128;
    int max_grid_size = 64;
    
    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 0;
    int max_iter = 100;
    int max_fmg_iter = 0;
    amrex::Real reltol = 1.e-11;

    amrex::Geometry geom;
    amrex::BoxArray grids;
    amrex::DistributionMapping dmap;

    amrex::MultiFab solution;
    amrex::MultiFab rhs;
    amrex::MultiFab exact_solution;
    amrex::MultiFab sigma;
};

#endif
