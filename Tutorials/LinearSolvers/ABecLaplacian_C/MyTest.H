#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MLMG.H>

class MyTest
{
public:

    MyTest ();

    void solve ();
    void writePlotfile () const;

private:

    void readParameters ();
    void initData ();
    void initProbPoisson ();
    void initProbABecLaplacian ();
    void solvePoisson ();
    void solveABecLaplacian ();

    int max_level = 1;
    int ref_ratio = 2;
    int n_cell = 128;
    int max_grid_size = 64;
    
    bool composite_solve = true;

    int prob_type = 1;  // 1. Poisson,  2. ABecLaplacian

    // For MLMG solver
    int verbose = 2;
    int cg_verbose = 0;
    int max_iter = 100;
    int max_fmg_iter = 0;
    int linop_maxorder = 2;
    bool agglomeration = true;
    bool consolidation = true;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    amrex::Vector<amrex::MultiFab> solution;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> exact_solution;
    amrex::Vector<amrex::MultiFab> acoef;
    amrex::Vector<amrex::MultiFab> bcoef;

    amrex::Real ascalar = 1.e-3;
    amrex::Real bscalar = 1.0;
};

#endif
