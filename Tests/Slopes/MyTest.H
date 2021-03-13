#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_Config.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_Array.H>

class MyTest
{
public:

    MyTest ();
    ~MyTest ();

    void compute_gradient ();
    void writePlotfile ();
    void initData ();

private:

    void initializeEB ();
    void readParameters ();
    void initGrids ();

    void initializeLinearData (int ilev);
    void initializeLinearDataFor2D (int ilev);
    void initializeLinearDataFor3D (int ilev);

    int max_level = 0;
    int ref_ratio = 2;
    int n_cell = 128;
    int max_grid_size = 64;
    amrex::Vector<int> is_periodic = amrex::Vector<int>(3,0);
    int eb_is_dirichlet = 0;
    int eb_is_homog_dirichlet = 0;
    amrex::Vector<amrex::Real> prob_lo = amrex::Vector<amrex::Real>(3,0.0);
    amrex::Vector<amrex::Real> prob_hi = amrex::Vector<amrex::Real>(3,1.0);

    std::string plot_file_name{"plot"};

    amrex::Vector<amrex::Real> scalars;

    int max_coarsening_level = 0;

    bool use_linear_1d = false;
    bool linear_1d_askew = false;
    amrex::Vector<amrex::Real> linear_1d_pt_on_top_wall = amrex::Vector<amrex::Real>(3,0.0);
    amrex::Real linear_1d_height = 1.0;
    amrex::Real linear_1d_rotation = 0.;
    int linear_1d_flow_dir = 0;
    int linear_1d_height_dir = 1;
    amrex::Real linear_1d_bottom = 0.0;
    amrex::Vector<amrex::Real> linear_1d_askew_rotation = amrex::Vector<amrex::Real>(2,0.0);
    int linear_1d_no_flow_dir = 0;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;
    amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory> > factory;

    amrex::Vector<amrex::MultiFab> phi;
    amrex::Vector<amrex::MultiFab> grad_x;
    amrex::Vector<amrex::MultiFab> grad_x_analytic;
    amrex::Vector<amrex::MultiFab> grad_y;
    amrex::Vector<amrex::MultiFab> grad_y_analytic;
    amrex::Vector<amrex::MultiFab> grad_z;
    amrex::Vector<amrex::MultiFab> grad_z_analytic;
    amrex::Vector<amrex::MultiFab> ccentr;

};

#endif
