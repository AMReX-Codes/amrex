#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLMG.H>
#include <AMReX_Array.H>

class MyTest
{
public:

    MyTest ();

    void apply ();
    void compute_gradient ();
    void writePlotfile ();
    void initData ();

private:

    void initializeEB ();
    void readParameters ();
    void initGrids ();

    void initializePoiseuilleData (int ilev);
    void initializePoiseuilleDataFor2D (int ilev);
    void initializePoiseuilleDataFor3D (int ilev);

    void initializeTriangleWaveData (int ilev);
    void initializeTriangleWaveDataFor2D (int ilev);

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

    amrex::Vector<amrex::BCRec> m_bcrec;

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
    bool use_poiseuille = false;
    bool use_triangle_wave = false;
    bool poiseuille_askew = false;
    amrex::Vector<amrex::Real> poiseuille_pt_on_top_wall = amrex::Vector<amrex::Real>(3,0.0);
    amrex::Real poiseuille_height = 1.0;
    amrex::Real poiseuille_rotation = 0.;
    int poiseuille_flow_dir = 0;
    int poiseuille_height_dir = 1;
    amrex::Real poiseuille_bottom = 0.0;
    amrex::Vector<amrex::Real> poiseuille_askew_rotation = amrex::Vector<amrex::Real>(2,0.0);
    int poiseuille_no_flow_dir = 0;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;
    amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory> > factory;

    amrex::Vector<amrex::MultiFab> phi;
    amrex::Vector<amrex::MultiFab> phi_ghost_resolved;
    amrex::Vector<amrex::MultiFab> phieb;
    amrex::Vector<amrex::MultiFab> grad_x;
    amrex::Vector<amrex::MultiFab> grad_x_analytic;
    amrex::Vector<amrex::MultiFab> grad_y;
    amrex::Vector<amrex::MultiFab> grad_y_analytic;
    amrex::Vector<amrex::MultiFab> grad_z;
    amrex::Vector<amrex::MultiFab> grad_z_analytic;
    amrex::Vector<amrex::MultiFab> grad_eb;
    amrex::Vector<amrex::MultiFab> grad_eb_analytic;
    amrex::Vector<amrex::MultiFab> ccentr;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> acoef;
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> > bcoef;
    amrex::Vector<amrex::MultiFab> bcoef_eb;
    amrex::Vector<amrex::MultiFab> lap_analytic;

    amrex::Vector<amrex::BCRec> const& get_bcrec () const noexcept { return m_bcrec; }

};

#endif
