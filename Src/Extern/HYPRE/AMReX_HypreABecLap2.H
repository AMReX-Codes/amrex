#ifndef AMREX_HYPREABECLAP2_H_
#define AMREX_HYPREABECLAP2_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MacBndry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>

#include "_hypre_sstruct_mv.h"
#include "_hypre_sstruct_ls.h"

namespace amrex
{

class HypreABecLap2
{
public:

    HypreABecLap2 (const BoxArray& grids,
                   const DistributionMapping& dmap,
                   const Geometry& geom_,
                   MPI_Comm comm_);
    ~HypreABecLap2 ();

    void setScalars (Real sa, Real sb);
    void setACoeffs (const MultiFab& alpha);
    void setBCoeffs (const Array<const MultiFab*,BL_SPACEDIM>& beta);
    void setVerbose (int _verbose);
    void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
                int max_iter, const BndryData& bndry, int max_bndry_order);

private :

    static constexpr int stencil_size = 2*AMREX_SPACEDIM+1;

    MPI_Comm comm = MPI_COMM_NULL;
    Geometry geom;

    int verbose = 0;

    int is_periodic[AMREX_SPACEDIM];

    MultiFab acoefs;
    Array<MultiFab,BL_SPACEDIM> bcoefs;
    Real scalar_a, scalar_b;

    HYPRE_SStructGrid     hgrid = NULL;
    HYPRE_SStructStencil  stencil = NULL;
    HYPRE_SStructGraph    graph = NULL;
    HYPRE_SStructMatrix   A = NULL;
    HYPRE_SStructVector   b = NULL;
    HYPRE_SStructVector   x = NULL;
    HYPRE_Solver          solver = NULL;

    Real rel_tol, abs_tol;
    int max_iter;

    void loadMatrix (const BndryData& bndry, int max_bndry_order);
    void finalizeMatrix ();
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
    void finalizeVectors ();

    void setupSolver (Real rel_tol_, Real abs_tol_, int max_iter_);
    void solveDoIt ();
    void clearSolver ();

    void getSolution (MultiFab& soln);
};

}

#endif
