#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLMG.H>
#include <AMReX_Array.H>

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
    int is_periodic = 0;
    int eb_is_dirichlet = 0;
    
    std::string plot_file_name{"plot"};

    amrex::Vector<amrex::Real> scalars;

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 2;
    int max_iter = 100;
    int max_fmg_iter = 0;
    int max_bottom_iter = 1000;
    amrex::Real reltol        = 1.e-12;
    amrex::Real bottom_reltol = 1.e-4;
    int linop_maxorder = 3;
    int max_coarsening_level = 30;
    bool use_hypre = false;
    bool use_petsc = false;
    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;
    amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory> > factory;

    amrex::Vector<amrex::MultiFab> phi;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> acoef;
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> > bcoef;
    amrex::Vector<amrex::MultiFab> bcoef_eb;
};

#endif
