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
    void writePlotfile () const;

    void initData ();

private:

    void readParameters ();

    int max_level = 0;
    int ref_ratio = 2;
    int n_cell = 128;
    int max_grid_size = 64;

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 0;
    int max_iter = 500;
    int max_fmg_iter = 0;
    int max_coarsening_level = 30;
    int use_hypre = 0;
    amrex::Real reltol = 1.e-11;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    amrex::Vector<amrex::MultiFab> solution;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> exact_solution;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> beta = {{AMREX_D_DECL(0.9,0.0,0.0)}};
};

#endif
