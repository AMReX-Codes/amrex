#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_Array.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MLMG.H>

// Solve the standard manufactured nodal Poisson problem with EB enabled.
// The MLEBNodeFDLaplacian default unit coefficient is used.

class MyTest
{
public:

    MyTest ();

    void solve ();
    void compute_norms () const;
    void writePlotfile () const;

    void initData ();

    [[nodiscard]] static int getNumTrials () { return 1; }
    [[nodiscard]] int getDoPlots () const { return do_plots; }

private:

    void readParameters ();

    int n_cell = 64;
    int max_grid_size = 32;
    int max_coarsening_level = 100;

    int verbose = 2;
    int bottom_verbose = 0;
    int max_iter = 200;
    int max_fmg_iter = 0;
    amrex::Real reltol = 1.e-10;

    int gpu_regtest = 0;
    bool do_plots = true;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    std::unique_ptr<amrex::EBFArrayBoxFactory> factory;

    amrex::Vector<amrex::MultiFab> solution;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> exact_solution;
    amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> grad_solution;
};

#endif
