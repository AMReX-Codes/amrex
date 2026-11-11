#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLMG.H>
#include <AMReX_Array.H>
#include <string>

class MyTest
{
public:

    MyTest ();

    void solve ();

    void writePlotfile ();

    void initData ();

    // public for cuda
    void compareWithMLHypre (amrex::MLEBABecLap& mleb, amrex::MLMG& mlmg,
                             amrex::Vector<amrex::MultiFab> const& bcdata,
                             amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& lobc,
                             amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> const& hibc);

private:

    void initializeEB ();
    void readParameters ();
    void initGrids ();
    void addFineGrids ();

    int max_level = 0;
    int ref_ratio = 2;
    int n_cell = 128;
    int max_grid_size = 64;

    int prob_type = 2;

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
    int agg_grid_size = -1;
    int con_grid_size = -1;

    int composite_solve = 0;
    bool eb_neumann = false;    // homogeneous Neumann on the EB instead of Dirichlet
    amrex::Real acoef_value = 0.0;
    amrex::Real ghost_value = 0.0; // in phi's ghost cells off the domain boundary
    amrex::Vector<int> is_periodic{AMREX_D_DECL(0,0,0)};
    // Region of level l: [fine_lo + (l-1)*fine_shrink, fine_hi - (l-1)*fine_shrink]
    amrex::Vector<amrex::Real> fine_lo{AMREX_D_DECL(0.,0.,0.)};
    amrex::Vector<amrex::Real> fine_hi{AMREX_D_DECL(1.,1.,1.)};
    amrex::Real fine_shrink = 0.0;
    // Removed from the finer levels (empty: nothing): an L-shaped fine level
    amrex::Vector<amrex::Real> fine_notch_lo;
    amrex::Vector<amrex::Real> fine_notch_hi;
    // b = 1 on the finer levels (MLMG gets it), given to HypreMLABecLap as
    // null (1: a null b means 1) or explicitly (2); the two must give the
    // same solution
    int null_fine_bcoef = 0;
    // HypreMLABecLap gets no a at all (null, that is zero; acoef must be 0)
    int null_acoef = 0;
    std::string domain_bc{"dirichlet"}; // or neumann, on the non-periodic domain boundaries
    amrex::Vector<amrex::Real> prob_hi{AMREX_D_DECL(1.,1.,1.)}; // the domain is [0, prob_hi]
    // Also solve with HypreMLABecLap (composite_solve): mlhypre_num_setups
    // setups with reuse and scaled coefficients and rhs (same solution), each
    // compared with MLMG's solution (all levels, mlhypre_max_diff) and
    // checked with MLMG's composite residual times the volume fraction
    // (finest level, and the coarser levels away from the coarse/fine
    // boundary, mlhypre_max_fine_res) and times its square root, as the
    // scaled rows the solver converges (the same cells,
    // mlhypre_max_res_scaled), with the EB average of the finer level's
    // solution (coarser levels, mlhypre_max_avg_diff) and against the exact
    // solution (error at most mlhypre_max_exact_ratio times MLMG's, and at
    // most mlhypre_max_exact_err, except for singular problems); negative
    // bounds are not checked
    bool use_mlhypre = false;
    int mlhypre_num_setups = 1;
    amrex::Real mlhypre_max_diff = -1.0;
    amrex::Real mlhypre_max_fine_res = -1.0;
    amrex::Real mlhypre_max_avg_diff = 1.e-12;
    amrex::Real mlhypre_max_exact_ratio = 10.0;
    amrex::Real mlhypre_max_res_scaled = -1.0;
    amrex::Real mlhypre_max_exact_err = -1.0;
    // NaN in the cells HypreMLABecLap does not use (EB-covered, covered by
    // a finer level) of its rhs and initial guess
    bool nan_in_covered = false;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;
    amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory> > factory;

    amrex::Vector<amrex::MultiFab> phi;
    amrex::Vector<amrex::MultiFab> phiexact;
    amrex::Vector<amrex::MultiFab> phieb;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> acoef;
    amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> > bcoef;
    amrex::Vector<amrex::MultiFab> bcoef_eb;
    amrex::Vector<amrex::MultiFab> fluxeb_sol;
    amrex::Vector<amrex::MultiFab> fluxeb_exact;
};

#endif
