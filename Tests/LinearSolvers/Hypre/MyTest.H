#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MultiFab.H>

class MyTest
{
public:

    MyTest ();
    void solve ();
    void initData ();

private:

    void initializeEB ();
    void readParameters ();
    void initGrids ();

    int n_cell = 128;
    int max_grid_size = 64;

    amrex::Real phi_domain = 0.0;
    amrex::Real phi_eb = 1.0;

    int verbose = 2;
    int max_iter = 100;
    amrex::Real reltol = 1.e-11;

    amrex::Geometry geom;
    amrex::BoxArray grids;
    amrex::DistributionMapping dmap;
    std::unique_ptr<amrex::EBFArrayBoxFactory> factory;

    amrex::MultiFab phi;
    amrex::MultiFab rhs;
};

#endif
