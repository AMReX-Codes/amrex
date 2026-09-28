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
    [[nodiscard]] amrex::Vector<std::string> const& getMultigridTypes () const {
        return multigrid_types;
    }
    void setMultigridType (std::string const& name);

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
#ifdef AMREX_USE_FLOAT
    amrex::Real reltol = 1.e-4F;
#else
    amrex::Real reltol = 1.e-11;
#endif

    bool agglomeration = true;
    bool consolidation = true;
    bool semicoarsening = false;
#if (AMREX_SPACEDIM == 2)
    bool rz = false;
#endif
    int max_coarsening_level = 30;
    int max_semicoarsening_level = 0;
    amrex::Vector<std::string> multigrid_types{"geometric"};
    amrex::MultigridType multigrid_type = amrex::MultigridType::geometric;

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
