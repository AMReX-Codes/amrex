#ifndef MY_TEST_H_
#define MY_TEST_H_

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif
#include <AMReX_MLMG.H>

#include <memory>
#include <string>

// Solve `div (sigma grad phi) = rhs` using MLEBNodeFDLaplacian.

class MyTest
{
public:

    MyTest ();

    void run ();
    void initData ();

private:

    void readParameters ();
    void solve (bool variable_sigma);
    void computeNorms (std::string const& label) const;

    int n_cell = 64;
    int max_grid_size = 32;

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 0;
    int max_iter = 100;
    int max_fmg_iter = 0;
    amrex::Real reltol = 1.e-11;

    bool agglomeration = true;
    bool consolidation = true;
    bool semicoarsening = false;
#if (AMREX_SPACEDIM == 2)
    bool rz = false;
#endif
    int max_coarsening_level = 30;
    int max_semicoarsening_level = 0;

    amrex::Geometry geom;
    amrex::BoxArray grids;
    amrex::DistributionMapping dmap;
#ifdef AMREX_USE_EB
    std::unique_ptr<amrex::EBFArrayBoxFactory> factory;
#endif

    amrex::MultiFab solution;
    amrex::MultiFab rhs;
    amrex::MultiFab exact_solution;
    amrex::MultiFab sigma;
};

#endif
